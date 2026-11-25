#' Distances from one tree to many, with the trees prepared once
#'
#' `DistanceOracle()` prepares a set of trees once, returning a function that
#' gives the distances from any one of them to all of them.
#'
#' Algorithms such as farthest-first sampling request one column of a distance
#' matrix at a time, each time comparing a single tree with the whole set.
#' Calling `Func(trees[[i]], trees)` re-reads the splits of every tree in each
#' call; the oracle reads them only once, so each column is cheaper.
#' `oracle(i)` returns the same values as `Func(trees[[i]], trees, ...)`, in the
#' shape that `Coreset::FarFirst(d = )` accepts.
#'
#' `NyeSimilarity()` returns a similarity unless `similarity = FALSE`.
#'
#' The prepared trees live in C++ memory, and are rebuilt if lost, as after
#' [`saveRDS()`]/[`readRDS()`].
#' Where trees do not share a set of four or more leaves, or would occupy more
#' than `maxBytes`, `oracle(i)` simply calls `Func(trees[[i]], trees, ...)`.
#' Columns are computed using `getOption("mc.cores")` threads.
#'
#' @param trees List of trees of class `phylo` or `Splits`, or a `multiPhylo`
#' object.
#' @param Func Function giving the distance between trees: one of
#' [`ClusteringInfoDistance()`], [`InfoRobinsonFoulds()`], [`RobinsonFoulds()`],
#' [`MatchingSplitDistance()`], [`NyeSimilarity()`] or
#' [`JaccardRobinsonFoulds()`].
#' @param \dots Further arguments to `Func`, such as `normalize`.
#' @param maxBytes Numeric giving the most memory, in bytes, that the prepared
#' trees may occupy.
#'
#' @return `DistanceOracle()` returns a function of `i`, an integer vector
#' indexing `trees`.  It returns the values of `Func` between `trees[[i]]` and
#' each tree in `trees`: a vector if `i` has length one, otherwise a matrix
#' with one column per entry of `i`.
#'
#' @examples
#' library("TreeTools", quietly = TRUE)
#' trees <- as.phylo(0:9, nTip = 12)
#' oracle <- DistanceOracle(trees, ClusteringInfoDistance)
#' oracle(3)
#' all.equal(oracle(3), ClusteringInfoDistance(trees[[3]], trees))
#' @template MRS
#' @family tree distances
#' @export
DistanceOracle <- function(trees, Func = ClusteringInfoDistance, ...,
                           maxBytes = 2^31) {
  if (inherits(trees, c("phylo", "Splits"))) {
    stop("`trees` must be a list of trees")
  }
  spec <- .OracleSpec(Func)
  args <- .OracleArgs(Func, ...)
  if (isTRUE(args[["reportMatching"]])) {
    stop("`reportMatching` is not supported")
  }
  normalize <- args[["normalize"]]
  if (spec[["name"]] == "MSD" &&
      (isTRUE(normalize) || is.function(normalize))) {
    stop("Please specify a function to generate a normalizing constant")
  }
  nTree <- length(trees)

  splits <- .OracleSplits(trees, maxBytes)
  if (is.null(splits)) {
    dots <- list(...)
    .Oracle(function(i) do.call(Func, c(list(trees[[i]], trees), dots)), nTree)
  } else {
    .PreparedOracle(splits, spec, args, names(trees))
  }
}

# Function of `i` that validates its input and applies `Column` to each entry
.Oracle <- function(Column, nTree) {
  force(Column)
  force(nTree)
  function(i) {
    if (!is.numeric(i) || length(i) == 0 || anyNA(i) ||
        any(i < 1 | i > nTree)) {
      stop("`i` must index trees 1 to ", nTree)
    }
    if (length(i) == 1) {
      Column(i)
    } else {
      vapply(i, Column, double(nTree))
    }
  }
}

# Splits of `trees` on a common leaf order, or NULL if they cannot be prepared
#' @importFrom TreeTools TipLabels
.OracleSplits <- function(trees, maxBytes) {
  labs <- TipLabels(trees)
  if (is.list(labs)) {
    if (!all(vapply(labs[-1], setequal, logical(1), labs[[1]]))) {
      return(NULL)
    }
    labs <- labs[[1]]
  }
  nTip <- length(labs)
  if (nTip < 4) {
    return(NULL)
  }
  .CheckMaxTips(nTip)
  splits <- .SplitsOn(trees, labs)
  nSplit <- vapply(splits, nrow, integer(1))
  # Rows and canonical forms, plus per-split and per-tree bookkeeping
  bytes <- sum(nSplit) * (16 * ceiling(nTip / 64) + 16) + 64 * length(trees)
  if (bytes > maxBytes) {
    return(NULL)
  }
  attr(splits, "nTip") <- nTip
  splits
}

.PreparedOracle <- function(splits, spec, args, treeNames) {
  nTip <- attr(splits, "nTip")
  Post <- .OraclePost(spec, args, nTip)
  k <- as.double(if (is.null(args[["k"]])) 1 else args[["k"]])
  allowConflict <- !isFALSE(args[["allowConflict"]])
  ptr <- NULL
  info <- NULL
  .Threads <- function() as.integer(getOption("mc.cores", 1L))
  .Prepare <- function() {
    ptr <<- cpp_oracle_prepare(splits, nTip, spec[["metric"]], spec[["info"]],
                               .Threads())
    info <<- attr(ptr, "info")
  }
  .Column <- function(i) {
    score <- cpp_oracle_column(ptr, i, k, allowConflict, .Threads())
    if (is.null(score)) {
      .Prepare()
      score <- cpp_oracle_column(ptr, i, k, allowConflict, .Threads())
    }
    names(score) <- treeNames
    Post(score, info[[i]], info)
  }
  .Prepare()
  .Oracle(.Column, length(splits))
}

# Kernel (a PrunedMetric) and tree information (a PrunedInfo) for `Func`
.OracleSpec <- function(Func) {
  if (identical(Func, ClusteringInfoDistance)) {
    list(name = "CID", metric = 1L, info = 1L)
  } else if (identical(Func, InfoRobinsonFoulds)) {
    list(name = "IRF", metric = 2L, info = 2L)
  } else if (identical(Func, RobinsonFoulds)) {
    list(name = "RF", metric = 7L, info = 3L)
  } else if (identical(Func, MatchingSplitDistance)) {
    list(name = "MSD", metric = 3L, info = 0L)
  } else if (identical(Func, NyeSimilarity)) {
    list(name = "Nye", metric = 6L, info = 3L)
  } else if (identical(Func, JaccardRobinsonFoulds)) {
    list(name = "JRF", metric = 6L, info = 3L)
  } else {
    stop("`Func` must be one of ClusteringInfoDistance, InfoRobinsonFoulds, ",
         "RobinsonFoulds, MatchingSplitDistance, NyeSimilarity or ",
         "JaccardRobinsonFoulds")
  }
}

# Values of each argument of `Func`, with defaults evaluated as `Func` would
.OracleArgs <- function(Func, ...) {
  ArgsOf <- Func
  body(ArgsOf) <- quote(
    mget(setdiff(names(formals(sys.function())), c("tree1", "tree2")))
  )
  ArgsOf(NULL, NULL, ...)
}

# Function converting kernel scores between tree `i`, with information `qi`,
# and all trees, with information `info`, to the values `Func` returns,
# following the arithmetic of `Func`'s one-to-many path on `nTip` leaves
.OraclePost <- function(spec, args, nTip) {
  normalize <- args[["normalize"]]
  similarity <- args[["similarity"]]
  switch(
    spec[["name"]],
    CID = {
      how <- if (is.function(normalize)) TRUE else normalize
      function(s, qi, info) {
        both <- qi + info
        ret <- .FloorNumericalNoise(both - s - s, both)
        .OracleNormalize(ret, qi, info, how, `+`)
      }
    },
    IRF = if (similarity) {
      function(s, qi, info) .OracleNormalize(s * 2, qi, info, normalize, `+`)
    } else {
      function(s, qi, info) {
        both <- qi + info
        ret <- .FloorNumericalNoise(both - s - s, both)
        .OracleNormalize(ret, qi, info, normalize, `+`)
      }
    },
    RF = function(s, qi, info) {
      ret <- if (similarity) (qi + info) - s else s
      .OracleNormalize(ret, qi, info, normalize, `+`)
    },
    MSD = function(s, qi, info) .OracleNormalize(s, qi, info, normalize, `+`),
    Nye = if (similarity) {
      normalizeMax <- args[["normalizeMax"]]
      function(s, qi, info) {
        if (normalizeMax) {
          qi <- nTip - 3
          info <- rep(qi, length(info))
        }
        .OracleNormalize(s, qi, info, normalize, .MeanOfTwo)
      }
    } else {
      function(s, qi, info) {
        .OracleNormalize((qi + info) - (s + s), qi, info, normalize, `+`)
      }
    },
    JRF = function(s, qi, info) {
      ret <- s * 2L
      if (!similarity) {
        ret <- (qi + info) - ret
      }
      .OracleNormalize(ret, qi, info, normalize, `+`)
    }
  )
}

# NormalizeInfo() for one tree, with information `qi`, against many, with
# information `info`; `Combine` must be vectorized
.OracleNormalize <- function(x, qi, info, how, Combine) {
  if (is.logical(how)) {
    if (how) x / Combine(qi, info) else x
  } else if (is.function(how)) {
    x / unlist(.mapply(how, list(qi, info), NULL))
  } else {
    x / how
  }
}
