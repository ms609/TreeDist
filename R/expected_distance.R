#' Expected information-theoretic tree distances under a random model
#'
#' Mean of the normalized [`ClusteringInfoDistance()`] (`Expected...Clustering...`)
#' or [`PhylogeneticInfoDistance()`] (`Expected...Phylogenetic...`) between a
#' random pair of `n`-leaf binary trees.  Used by
#' [`AdjustedClusteringInfoDistance()`] and
#' [`AdjustedPhylogeneticInfoDistance()`] to chance-adjust their scores
#' following @SmithDist.
#'
#' ## Size-conditioned null (`null = "size"`, default)
#'
#' Three back-ends are available:
#'
#' \describe{
#'   \item{`"lookup"`}{Reads the mean from the precomputed reference table
#'     [`randomTreeDistances`], shipped with \pkg{TreeDist} for `n` in
#'     `4:200`.  If `n` falls outside this range, the call silently falls
#'     back to `"mc"` and a message is emitted.}
#'   \item{`"enumerate"`}{Exact mean over every unordered pair of binary
#'     trees on `n` tips.  Capped at `n <= 7` because the
#'     `((2n - 3)!! choose 2)` pair count grows super-exponentially
#'     (135 135 trees for `n = 8`).}
#'   \item{`"mc"`}{Monte Carlo estimate over `nSim` random pairs drawn with
#'     `TreeTools::RandomTree(n, root = TRUE)`, matching the protocol used
#'     to build [`randomTreeDistances`].}
#' }
#'
#' The functions vectorise over `n`: passing a vector returns a vector of
#' expected distances.
#'
#' ## Shape-conditioned null (`null = "shape"`)
#'
#' Conditions on the unlabelled shape of both trees, so that only tip labels
#' are randomised.  Two trees with the same shapes as `tree1` and `tree2` are
#' compared, but with labels drawn independently at random.  This null is
#' strictly stronger than the size null: it removes variance attributable to
#' tree shape, leaving only variance from label placement.
#'
#' Two back-ends are available when `null = "shape"`:
#'
#' \describe{
#'   \item{`"shape-lookup"`}{Looks up `(shape1, shape2)` in the precomputed
#'     table [`shapeExpectedDistances`] (shipped for `n` in `4:10`; falls
#'     back to `"shape-mc"` with a message for larger `n`).}
#'   \item{`"shape-mc"`}{Label-shuffle Monte Carlo: independently samples
#'     `nSim` random label assignments for each shape via
#'     `TreeTools::RootedTreeWithShape()` and returns the mean distance.  The
#'     standard deviation is attached as attribute `"sd"`.}
#' }
#'
#' When `null = "shape"`, both `tree1` and `tree2` must be supplied as single
#' `phylo` objects (not lists).
#'
#' @param n Integer; number of tips.  Must be at least 4.  May be a vector.
#'   Ignored when `null = "shape"` (tip count is taken from `tree1`).
#' @param method Character; back-end to use.  One of `"lookup"`,
#'   `"enumerate"`, `"mc"` (for `null = "size"`), or `"shape-lookup"`,
#'   `"shape-mc"` (for `null = "shape"`).  When `null = "shape"` and
#'   `method` is not one of the shape-specific values, it defaults to
#'   `"shape-lookup"`.
#' @param nSim Integer; Monte Carlo replicates when `method = "mc"`,
#'   `"shape-mc"`, or the lookup falls back to MC.
#' @param null Character; null model.  `"size"` (default) fixes only the
#'   number of tips; `"shape"` additionally fixes the unlabelled tree shape.
#' @param tree1,tree2 Single `phylo` objects; required when `null = "shape"`.
#'   Ignored when `null = "size"`.
#' @param \dots Currently unused; reserved for future tuning parameters.
#'
#' @returns A numeric scalar (length-1 vector) of the normalized expected
#'   distance.  When `null = "size"`, the function vectorises over `n`.
#'   When `null = "shape"` and `method = "shape-mc"`, the return value
#'   carries an attribute `"sd"` giving the Monte Carlo standard deviation.
#'
#' @examples
#' # Fast lookup (shipped table)
#' ExpectedClusteringInfoDistance(8)
#' ExpectedPhylogeneticInfoDistance(c(8, 16, 32))
#'
#' # Explicit enumeration for very small n
#' ExpectedClusteringInfoDistance(6, method = "enumerate")
#'
#' # On-the-fly Monte Carlo with a tiny replicate budget
#' ExpectedPhylogeneticInfoDistance(8, method = "mc", nSim = 50L)
#'
#' # Shape-conditioned null (lookup)
#' library("TreeTools", quietly = TRUE)
#' t1 <- PectinateTree(8); t2 <- BalancedTree(8)
#' ExpectedClusteringInfoDistance(null = "shape", tree1 = t1, tree2 = t2)
#'
#' # Shape-conditioned null (Monte Carlo)
#' ExpectedClusteringInfoDistance(null = "shape", method = "shape-mc",
#'                                tree1 = t1, tree2 = t2, nSim = 200L)
#'
#' @template MRS
#' @references
#' \insertAllCited{}
#' @encoding UTF-8
#' @family tree distances
#' @seealso [`randomTreeDistances`], [`shapeExpectedDistances`]
#' @importFrom TreeTools NRootedShapes NTip Preorder RootedTreeShape
#'   RootedTreeWithShape RandomTree
#' @export
ExpectedClusteringInfoDistance <- function(n,
                                           method = c("lookup", "enumerate",
                                                      "mc", "shape-lookup",
                                                      "shape-mc"),
                                           nSim = 1e4L,
                                           null = c("size", "shape"),
                                           tree1 = NULL, tree2 = NULL,
                                           ...) {
  null <- match.arg(null)
  if (null == "shape") {
    .ExpectedShapeDistance(tree1 = tree1, tree2 = tree2,
                           method = method, nSim = nSim,
                           dist_fun = ClusteringInfoDistance,
                           table_col_mean = "mean_cid",
                           table_col_sd   = "sd_cid")
  } else {
    .ExpectedRandomDistance(n, method = method, nSim = nSim,
                            dist_fun = ClusteringInfoDistance,
                            table_col = "cid_mean")
  }
}

#' @rdname ExpectedClusteringInfoDistance
#' @export
ExpectedPhylogeneticInfoDistance <- function(n,
                                             method = c("lookup", "enumerate",
                                                        "mc", "shape-lookup",
                                                        "shape-mc"),
                                             nSim = 1e4L,
                                             null = c("size", "shape"),
                                             tree1 = NULL, tree2 = NULL,
                                             ...) {
  null <- match.arg(null)
  if (null == "shape") {
    .ExpectedShapeDistance(tree1 = tree1, tree2 = tree2,
                           method = method, nSim = nSim,
                           dist_fun = PhylogeneticInfoDistance,
                           table_col_mean = "mean_pid",
                           table_col_sd   = "sd_pid")
  } else {
    .ExpectedRandomDistance(n, method = method, nSim = nSim,
                            dist_fun = PhylogeneticInfoDistance,
                            table_col = "pid_mean")
  }
}

# Internal: shared back-end for the Expected... functions (size-conditioned).
.ExpectedRandomDistance <- function(n, method, nSim, dist_fun, table_col) {
  method <- match.arg(method, c("lookup", "enumerate", "mc",
                                "shape-lookup", "shape-mc"))
  # Ignore shape-specific methods if accidentally passed here
  if (method %in% c("shape-lookup", "shape-mc")) method <- "lookup"
  n <- as.integer(n)
  if (anyNA(n) || any(n < 4L)) {
    stop("`n` must be an integer vector with all values >= 4.")
  }
  nSim <- as.integer(nSim)
  if (length(nSim) != 1L || is.na(nSim) || nSim < 2L) {
    stop("`nSim` must be a positive integer (>= 2).")
  }

  vapply(n, function(ni) {
    switch(method,
           lookup = .ExpectedLookup(ni, nSim = nSim, dist_fun = dist_fun,
                                    table_col = table_col),
           enumerate = .ExpectedEnumerate(ni, dist_fun = dist_fun),
           mc = .ExpectedMC(ni, nSim = nSim, dist_fun = dist_fun))
  }, numeric(1))
}

.ExpectedLookup <- function(n, nSim, dist_fun, table_col) {
  tbl <- get("randomTreeDistances", envir = asNamespace("TreeDist"))
  hit <- tbl[tbl[["n"]] == n, , drop = FALSE]
  if (nrow(hit) == 1L) {
    return(hit[[table_col]])
  }
  message("`n = ", n,
          "` is outside the shipped lookup range; using Monte Carlo with ",
          nSim, " replicates.")
  .ExpectedMC(n, nSim = nSim, dist_fun = dist_fun)
}

# Cap on enumeration. `n = 7` gives 10 395 rooted binary trees and
# ~5.4e7 unordered pairs; runs in minutes. `n = 8` (135 135 trees) is
# memory-prohibitive.
.enumerateCap <- 7L

.ExpectedEnumerate <- function(n, dist_fun) {
  if (n > .enumerateCap) {
    stop("`method = \"enumerate\"` is capped at n <= ", .enumerateCap,
         "; got n = ", n,
         ".  Use `method = \"mc\"` or `method = \"lookup\"` for larger n.")
  }
  nTrees <- TreeTools::NRooted(n)
  # S3 dispatch: `as.phylo.numeric` is registered by TreeTools (which is
  # imported), reached via the `ape` generic.
  trees <- ape::as.phylo(seq_len(nTrees) - 1L, nTip = n)
  # `dist_fun(list, normalize = TRUE)` returns a `dist` over the
  # `nTrees * (nTrees - 1) / 2` unordered off-diagonal pairs.  The
  # independent-draw Monte Carlo protocol that built `randomTreeDistances`
  # allows self-pairs (probability 1/nTrees), so the matching exact
  # expectation is over all `nTrees^2` ordered pairs *including* the
  # zero-distance diagonal:
  #   E_inc = (2 * sum(d_off)) / nTrees^2.
  (2 * sum(dist_fun(trees, normalize = TRUE))) / (nTrees * nTrees)
}

.ExpectedMC <- function(n, nSim, dist_fun) {
  dists <- vapply(seq_len(nSim), function(i) {
    t1 <- TreeTools::RandomTree(n, root = TRUE)
    t2 <- TreeTools::RandomTree(n, root = TRUE)
    dist_fun(t1, t2, normalize = TRUE)
  }, numeric(1))
  mean(dists)
}

# Internal: shape-conditioned expected distance dispatcher.
#
# Validates inputs then dispatches to shape-lookup or shape-mc back-end.
.ExpectedShapeDistance <- function(tree1, tree2, method, nSim,
                                   dist_fun, table_col_mean, table_col_sd) {
  if (is.null(tree1) || is.null(tree2)) {
    stop("Both `tree1` and `tree2` must be supplied when `null = \"shape\"`.")
  }
  if (!inherits(tree1, "phylo") || !inherits(tree2, "phylo")) {
    stop("`tree1` and `tree2` must be single `phylo` objects ",
         "when `null = \"shape\"`.")
  }

  # Normalise method: default to shape-lookup; treat size methods as
  # shape-lookup so that Adjusted* can pass method through transparently.
  method <- match.arg(method, c("lookup", "enumerate", "mc",
                                "shape-lookup", "shape-mc"))
  shape_method <- if (method == "shape-mc") "shape-mc" else "shape-lookup"

  n <- TreeTools::NTip(tree1)
  if (TreeTools::NTip(tree2) != n) {
    stop("`tree1` and `tree2` must have the same number of tips for ",
         "`null = \"shape\"`.")
  }

  s1 <- as.integer(TreeTools::RootedTreeShape(tree1))
  s2 <- as.integer(TreeTools::RootedTreeShape(tree2))
  lo <- min(s1, s2)
  hi <- max(s1, s2)

  nSim <- as.integer(nSim)
  if (length(nSim) != 1L || is.na(nSim) || nSim < 2L) {
    stop("`nSim` must be a positive integer (>= 2).")
  }

  if (shape_method == "shape-lookup") {
    .ExpectedShapeLookup(n = n, lo = lo, hi = hi, nSim = nSim,
                         dist_fun = dist_fun,
                         table_col_mean = table_col_mean,
                         table_col_sd   = table_col_sd)
  } else {
    .ExpectedShapeMC(n = n, s1 = s1, s2 = s2, nSim = nSim,
                     dist_fun = dist_fun)
  }
}

# Internal: look up (lo, hi) in shapeExpectedDistances for this n.
# Falls back to shape-mc with a message if n is outside the shipped range.
.ExpectedShapeLookup <- function(n, lo, hi, nSim, dist_fun,
                                 table_col_mean, table_col_sd) {
  tbl_list <- .ShapeExpectedDistances()
  key <- as.character(n)
  if (!key %in% names(tbl_list)) {
    message("`n = ", n, "` is outside the shipped shape lookup range ",
            "(n = 4:10); using shape Monte Carlo with ", nSim, " replicates.")
    return(.ExpectedShapeMC(n = n, s1 = lo, s2 = hi,
                            nSim = nSim, dist_fun = dist_fun))
  }
  tbl <- tbl_list[[key]]
  row <- tbl[tbl[["shape1"]] == lo & tbl[["shape2"]] == hi, , drop = FALSE]
  if (nrow(row) != 1L) {
    message("Shape pair (", lo, ", ", hi, ") not found in shipped table for ",
            "n = ", n, "; using shape Monte Carlo with ", nSim, " replicates.")
    return(.ExpectedShapeMC(n = n, s1 = lo, s2 = hi,
                            nSim = nSim, dist_fun = dist_fun))
  }
  row[[table_col_mean]]
}

# Internal: label-shuffle MC for shape-conditioned expected distance.
# Independently randomises tip labels for each shape; returns mean with
# attribute "sd".
# Note: RootedTreeWithShape() returns trees without an 'order' attribute;
# Preorder() is required before passing to as.Splits / dist functions.
.ExpectedShapeMC <- function(n, s1, s2, nSim, dist_fun) {
  lab_pool <- as.character(seq_len(n))
  dists <- vapply(seq_len(nSim), function(i) {
    t1 <- TreeTools::Preorder(
            TreeTools::RootedTreeWithShape(s1, n, sample(lab_pool)))
    t2 <- TreeTools::Preorder(
            TreeTools::RootedTreeWithShape(s2, n, sample(lab_pool)))
    dist_fun(t1, t2, normalize = TRUE)
  }, numeric(1))
  result <- mean(dists)
  attr(result, "sd") <- sd(dists)
  result
}

# Shared chance-adjustment engine used by ACID and APhID.
#
# Computes the normalized distance, then maps it to `1 - D / E`.
# `E` is either supplied via `expected` (scalar or matrix matching the shape
# of `D`) or looked up via `expected_fun` on a per-pair basis.  Mirrors the
# return shape (numeric / dist / matrix) of `dist_fun()`.
#
# When `null = "shape"`, single-pair dispatch is used (shape null only
# supports scalar tree inputs, not lists).
.AdjustInfoDistance <- function(tree1, tree2, expected, dist_fun,
                                expected_fun, null = "size", ...) {
  null <- match.arg(null, c("size", "shape"))

  d <- dist_fun(tree1, tree2, normalize = TRUE)
  d_attr <- attributes(d)

  e <- if (!is.null(expected)) {
    expected
  } else if (null == "shape") {
    # Shape null: single pair only.
    if (is.null(tree2)) {
      stop("`null = \"shape\"` requires `tree2` to be supplied.")
    }
    expected_fun(null = "shape", tree1 = tree1, tree2 = tree2, ...)
  } else {
    .ExpectedFor(tree1, tree2, expected_fun, ...)
  }

  # Chance-adjusted distance (D/E[D|n] convention):
  #   identity (D = 0)               -> 0
  #   random pair (D ~= E)           -> ~1
  #   more different than random     -> >1  (informative tail, not a bug)
  adj <- d / e
  degenerate <- e == 0
  if (any(degenerate)) {
    warning("Expected distance == 0 for some pair; returning NA there.")
    if (length(degenerate) == 1L) {
      adj[] <- NA_real_
    } else {
      adj[degenerate] <- NA_real_
    }
  }
  # Reattach attributes of the original distance object (preserves
  # `dist` class, matrix dims, etc.).
  attributes(adj) <- d_attr
  adj
}

# Build `E[D | n]` matched to the shape of `dist_fun(tree1, tree2)`.
# Returns a scalar when all relevant trees share the same `n`; otherwise a
# matrix or `dist` object of per-pair expectations, where the effective `n`
# for a pair is `min(NTip(t_i), NTip(t_j))` (CID/PID drops unshared tips
# before comparison).
.ExpectedFor <- function(tree1, tree2, expected_fun, ...) {
  n1 <- TreeTools::NTip(tree1)
  if (is.null(tree2)) {
    if (length(n1) == 1L || all(n1 == n1[[1L]])) {
      return(expected_fun(n1[[1L]], ...))
    }
    n_mat <- outer(n1, n1, pmin)
    e_mat <- matrix(expected_fun(as.integer(n_mat), ...),
                    nrow = nrow(n_mat), ncol = ncol(n_mat))
    return(stats::as.dist(e_mat))
  }
  n2 <- TreeTools::NTip(tree2)
  if (length(n1) == 1L && length(n2) == 1L && n1 == n2) {
    return(expected_fun(n1, ...))
  }
  if (all(n1 == n1[[1L]]) && all(n2 == n2[[1L]]) && n1[[1L]] == n2[[1L]]) {
    return(expected_fun(n1[[1L]], ...))
  }
  n_mat <- outer(n1, n2, pmin)
  matrix(expected_fun(as.integer(n_mat), ...),
         nrow = nrow(n_mat), ncol = ncol(n_mat))
}

# -------------------------------------------------------------------------
# Lazy-loaded shape lookup table
# -------------------------------------------------------------------------

# Cache environment for shapeExpectedDistances.
.shape_cache <- new.env(parent = emptyenv())

#' Load the shape-conditioned expected-distance lookup table
#'
#' Reads `inst/extdata/shapeExpectedDistances.rds` on first call and caches
#' the result.  Subsequent calls return the cached object.  The table covers
#' `n = 4:10`; for `n > 10` use `method = "shape-mc"`.
#'
#' @returns A named list, one element per `n` (named `"4"`, `"5"`, ...,
#'   `"10"`).  Each element is a data frame with columns `shape1`, `shape2`,
#'   `mean_cid`, `sd_cid`, `mean_pid`, `sd_pid`, `n_pairs`.
#'
#' @keywords internal
.ShapeExpectedDistances <- function() {
  if (!exists("tbl", envir = .shape_cache, inherits = FALSE)) {
    path <- system.file("extdata", "shapeExpectedDistances.rds",
                        package = "TreeDist")
    if (!nzchar(path)) {
      stop("shapeExpectedDistances.rds not found in TreeDist/inst/extdata. ",
           "Run data-raw/build_shape_lookup.R to build it.")
    }
    assign("tbl", readRDS(path), envir = .shape_cache)
  }
  get("tbl", envir = .shape_cache, inherits = FALSE)
}
