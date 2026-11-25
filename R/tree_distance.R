#' Generalized Robinson&ndash;Foulds distance
#' 
#' An internal function to calculate Generalized Robinson&ndash;Foulds
#' distances from splits.
#'
#' Note that no checks will be made to confirm that `splits1` and `splits2`
#' contain the same leaves in the same order.
#' This is the responsibility of the calling function.
#' 
#' @inheritParams SharedPhylogeneticInfoSplits
#' @param nTip Integer specifying the number of leaves in each split.
#' @param PairScorer function taking four arguments, `splits1`, `splits2`,
#' `nSplits1`, `nSplits2`, which should return the score of each pair of splits
#' in a two-dimensional matrix.  Additional parameters may be specified via 
#' \code{...}.
#' @param maximize Logical specifying whether the optimal matching maximizes
#' or minimizes the scores obtained by `PairScorer()`.
#' @param \dots Additional parameters to `PairScorer()`.
#' 
#' @return A numeric value specifying the score of the tree pairs under the 
#' specified pair scorer. If `reportMatching = TRUE`, attribute also list:
#' 
#' - `matching`: which split in `splits2` is optimally matched to each split in 
#'   `split1` (`NA` if not matched);
#' 
#' - `matchedSplits`: Textual representation of each match
#'
#' - `matchedScores`: Scores for matched split.
#' 
#' - `pairScores`: Calculated scores for each possible matching of each split.
#'
#' @keywords internal
#' @template MRS
#' @encoding UTF-8
#' @export
GeneralizedRF <- function(splits1, splits2, nTip, PairScorer, 
                           maximize, reportMatching, ...) {
  .ValidateSplitArgs(splits1, splits2, nTip)
  nSplits1 <- dim(splits1)[[1]]
  nSplits2 <- dim(splits2)[[1]]
  
  solution <- PairScorer(splits1, splits2, nTip, ...)
  ret <- solution[["score"]]
  
  if (reportMatching) {
    matching <- solution[["matching"]]
    matching[matching > nSplits2 | matching == 0L] <- NA
    if (nSplits1 < nSplits2) {
      matching <- matching[seq_len(nSplits1)]
    }
    attr(ret, "matching") <- matching
    
    # If reporting matching, we're not worried about performance
    pairScores <- matrix(0, nSplits1, nSplits2)
    for (i in seq_len(nSplits1)) {
      for (j in seq_len(nSplits2)) {
        pairScores[i, j] <- PairScorer(splits1[i, , drop = FALSE], 
                                       splits2[j, , drop = FALSE],
                                       nTip = nTip, ...)[["score"]]
      }
    }
    
    if (!is.null(attr(splits1, "tip.label"))) {
      matched1 <- !is.na(matching)
      matched2 <- matching[matched1]
      matched1 <- which(matched1)
      
      # Textual representation of matchings
      attr(ret, "matchedSplits") <- 
        ReportMatching(splits1[[matched1]],
                       splits2[[matched2]],
                       realMatch = if (maximize) {
                         pairScores[matrix(c(matched1, matched2), ncol = 2L)] > 0
                       } else rep(TRUE, length(matched1)))
    }
    
    attr(ret, "matchedScores") <- pairScores[
      matrix(c(seq_along(matching), matching), ncol = 2L)]
    
    attr(ret, "pairScores") <- pairScores
  }
  
  # Return:
  ret
}


.MaxValue <- function(tree1, tree2, Value) {
  lab1 <- TipLabels(tree1)
  sameTips <- .AllTipsSame(lab1, TipLabels(tree2))
  
  if (!sameTips) {
    # Trees in `tree1` do not all share leaves if !sameTips and tree2 is NULL
    if (!is.na(.PrunedInfo(Value))) {
      pruned <- .PrunedPairs(tree1, tree2, Value = Value)
      ret <- pruned[["info1"]] + pruned[["info2"]]
      single <- c("phylo", "Splits")
      if (!is.null(tree2) && !inherits(tree1, single) &&
          !inherits(tree2, single)) {
        dim(ret) <- c(length(tree1), length(tree2))
      }
      return(ret)
    }
    trees <- .SharedOnly(tree1, tree2)
    return(Value(trees[[1]]) + Value(trees[[2]]))
  }
  
  value1 <- Value(tree1)
  if (is.null(tree2)) {
    .PairwiseSums(value1)
  } else {
    outer(value1, Value(tree2), "+")[, , drop = TRUE]
  }
}

# Fast path for *Distance() functions: computes pairwise distances and
# per-tree info in a single as.Splits() conversion.  Returns NULL when the
# fast path is not applicable.
#' @importFrom TreeTools as.Splits TipLabels
.FastDistPath <- function(tree1, tree2, reportMatching,
                          cpp_batch_fn, cpp_entropy_fn) {
  if (!is.null(tree2) || reportMatching) return(NULL)
  if (inherits(tree1, c("phylo", "Splits"))) return(NULL)
  if (!is.null(getOption("TreeDist-cluster"))) return(NULL)
  
  labs <- TipLabels(tree1)
  if (is.list(labs)) {
    if (!all(vapply(labs[-1], setequal, logical(1), labs[[1]]))) return(NULL)
    tipLabels <- labs[[1]]
  } else {
    tipLabels <- labs
  }
  nTip <- length(tipLabels)
  if (nTip < 4) return(NULL)
  .CheckMaxTips(nTip)
  
  splits_list <- .SplitsOn(tree1, tipLabels)
  n_threads <- as.integer(getOption("mc.cores", 1L))
  
  info_vec <- cpp_batch_fn(splits_list, as.integer(nTip), n_threads)
  entropies <- cpp_entropy_fn(splits_list, as.integer(nTip))
  
  list(
    info = structure(info_vec, class = "dist",
                     Size = length(tree1), Labels = names(tree1),
                     Diag = FALSE, Upper = FALSE),
    entropies = entropies
  )
}

# Fast path for cross-pairs (ManyMany): avoids duplicate as.Splits() calls by
# computing both pairwise distances and per-tree entropies in a single pass.
# Returns NULL when not applicable. When applicable, returns list with:
#   $dists: nA × nB matrix of pairwise distances, or a vector if either
#           input is a single tree
#   $info1, $info2: each tree's information, if `withInfo`
#   $infoBoth: summed per-tree information, with the same shape as $dists
#' @importFrom stats setNames
#' @importFrom TreeTools as.Splits TipLabels
.FastManyManyPath <- function(tree1, tree2, reportMatching,
                              cpp_cross_pairs_fn, withInfo = TRUE) {
  if (is.null(tree2) || reportMatching) return(NULL)
  if (!is.null(getOption("TreeDist-cluster"))) return(NULL)
  x <- .CrossSplits(tree1, tree2)
  if (is.null(x)) return(NULL)
  n_threads <- as.integer(getOption("mc.cores", 1L))
  dists <- cpp_cross_pairs_fn(x[["splits1"]], x[["splits2"]], x[["nTip"]],
                              n_threads, with_info = withInfo)
  info1 <- attr(dists, "info_a")
  info2 <- attr(dists, "info_b")
  attributes(dists) <- list(dim = dim(dists))
  infoBoth <- if (withInfo) outer(info1, info2, "+")
  if (x[["single1"]]) {
    dists <- setNames(dists[1, ], x[["names2"]])
    infoBoth <- infoBoth[1, ]
  } else if (x[["single2"]]) {
    dists <- setNames(dists[, 1], x[["names1"]])
    infoBoth <- infoBoth[, 1]
  } else {
    dimnames(dists) <- list(x[["names1"]], x[["names2"]])
  }

  list(dists = dists, info1 = info1, info2 = info2, infoBoth = infoBoth)
}

# Splits for a cross-pairs kernel to compare `tree1` with `tree2`, at least
# one of which is a list; NULL unless all trees share a set of >= 4 leaves.
.CrossSplits <- function(tree1, tree2) {
  single1 <- inherits(tree1, c("phylo", "Splits"))
  single2 <- inherits(tree2, c("phylo", "Splits"))
  if (single1 && single2) return(NULL)
  names1 <- names(tree1)
  names2 <- names(tree2)
  if (single1) tree1 <- .WrapSingle(tree1)
  if (single2) tree2 <- .WrapSingle(tree2)

  if (single2) {
    lab2 <- TipLabels(tree2)
    lab1 <- .TipLabelsLike(tree1, lab2)
  } else {
    lab1 <- TipLabels(tree1)
    lab2 <- .TipLabelsLike(tree2, lab1)
  }
  
  # Check tip label agreement
  if (is.list(lab1)) {
    if (!all(vapply(lab1[-1], setequal, logical(1), lab1[[1]]))) return(NULL) # nocov
    tipLabels1 <- lab1[[1]] # nocov
  } else {
    tipLabels1 <- lab1
  }
  
  if (is.list(lab2)) {
    if (!all(vapply(lab2[-1], setequal, logical(1), lab2[[1]]))) return(NULL) # nocov
    tipLabels2 <- lab2[[1]] # nocov
  } else {
    tipLabels2 <- lab2
  }
  
  # Only use fast path if both collections have the same tip set
  if (!setequal(tipLabels1, tipLabels2)) return(NULL)
  
  nTip <- length(tipLabels1)
  if (nTip < 4) return(NULL)
  .CheckMaxTips(nTip)
  
  list(
    splits1 = .SplitsOn(tree1, tipLabels1),
    splits2 = .SplitsOn(tree2, tipLabels1),
    nTip = nTip, single1 = single1, single2 = single2,
    names1 = names1, names2 = names2
  )
}

# Information against which NormalizeInfo() should scale scores between trees
# carrying information `info1` and `info2`; NULL where `how` supplies it, or
# where no normalization is required.
.InfoToNormalize <- function(info1, info2, how, Combine, lower = FALSE) {
  if (!isTRUE(how) && !is.function(how)) {
    return(NULL)
  }
  single <- length(info1) == 1 || length(info2) == 1
  if (is.function(how) && single) {
    # User functions need not be vectorized, as in NormalizeInfo()
    return(unlist(.mapply(how, list(info1, info2), NULL)))
  }
  ret <- outer(info1, info2, if (is.function(how)) how else Combine)
  if (lower) {
    ret[lower.tri(ret)]
  } else if (single) {
    as.vector(ret)
  } else {
    ret
  }
}

# Similarities, and the information to normalize them against, from the batch
# kernels that score them; NULL where these do not apply.
.FastSimilarity <- function(tree1, tree2, normalize, reportMatching, diag,
                            cpp_all_pairs_fn, cpp_info_fn, cpp_cross_pairs_fn) {
  if (is.null(tree2)) {
    fast <- .FastDistPath(tree1, NULL, reportMatching, cpp_all_pairs_fn,
                          cpp_info_fn)
    if (is.null(fast)) {
      return(NULL)
    }
    info <- fast[["entropies"]]
    scores <- fast[["info"]]
    if (diag) {
      scores <- as.matrix(scores)
      diag(scores) <- info
    }
    infoInBoth <- .InfoToNormalize(info, info, normalize, .PairMean,
                                   lower = !diag)
  } else {
    fast <- .FastManyManyPath(tree1, tree2, reportMatching,
                              cpp_cross_pairs_fn,
                              isTRUE(normalize) || is.function(normalize))
    if (is.null(fast)) {
      return(NULL)
    }
    scores <- fast[["dists"]]
    infoInBoth <- .InfoToNormalize(fast[["info1"]], fast[["info2"]], normalize,
                                   .PairMean)
  }
  list(scores = scores, infoInBoth = infoInBoth)
}

# Lower-tri pairwise sums: outer(x, x, "+")[lower.tri(.)]
.PairwiseSums <- function(x) {
  g <- outer(x, x, "+")
  g[lower.tri(g)]
}

# Floor sub-noise distances to zero before normalization.
# Two sources of numerical noise scale with treesIndependentInfo:
#  (1) LAP int64 cost-matrix quantization in *Splits scoring; per-cell
#      truncation of up to (max_possible / BIG) bits, summed over n_splits.
#  (2) Float-accumulation drift between independently-built tables (e.g.
#      InfoRobinsonFoulds vs cpp_splitwise_info_batch sum the same per-split
#      info contributions, but using different lookup-table constructions).
# Both grow with the magnitude of the answer, so an absolute sqrt(eps)
# tolerance becomes too tight beyond a few thousand tips. Scaling by
# treesIndependentInfo self-adjusts; pmax(1, ·) preserves the original
# tolerance for tiny trees where these errors are negligible anyway.
.FloorNumericalNoise <- function(ret, treesIndependentInfo) {
  tol <- pmax(1, treesIndependentInfo) * .Machine[["double.eps"]] ^ 0.5
  ret[ret < tol] <- 0
  ret
}

.AllTipsSame <- function(x, y) {
  if (is.list(x)) {
    xPrime <- x[[1]]
    if (length(x) > 1 && !all(vapply(x[-1], setequal, logical(1), xPrime))) {
      return(FALSE)
    }
  } else {
    xPrime <- x
  }
  if (is.null(y)) {
    TRUE
  } else {
    if (is.list(y)) {
      all(vapply(y, setequal, logical(1), xPrime))
    } else {
      setequal(xPrime, y)
    }
  }
}

#' Are splits compatible?
#' 
#' Determine whether splits are compatible (concave); i.e. they can both occur
#' on a single tree.
#' 
#' @template split12Params
#' @return `SplitsCompatible()` returns a logical specifying whether the splits
#' provided are compatible with one another.
#' 
#' @examples 
#' A <- TRUE
#' B <- FALSE
#' SplitsCompatible(c(A, A, A, B, B, B),
#'                  c(A, A, B, B, B, B))
#' SplitsCompatible(c(A, A, A, B, B, B),
#'                  c(A, A, B, B, B, A))
#' @template MRS
#' @export
SplitsCompatible <- function(split1, split2) {
  # Return:
  (
    all (split1[split2]) ||
    all(split1[!split2]) ||
    all(!split1[split2]) ||
    all(!split1[!split2])
  )
}
