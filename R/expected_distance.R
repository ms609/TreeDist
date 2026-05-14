#' Expected information-theoretic tree distances under a random model
#'
#' Mean of the normalized [`ClusteringInfoDistance()`] (`Expected...Clustering...`)
#' or [`PhylogeneticInfoDistance()`] (`Expected...Phylogenetic...`) between a
#' random pair of `n`-leaf binary trees.  Used by
#' [`AdjustedClusteringInfoDistance()`] and
#' [`AdjustedPhylogeneticInfoDistance()`] to chance-adjust their scores
#' following @SmithDist.
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
#' @param n Integer; number of tips.  Must be at least 4.  May be a vector.
#' @param method Character; back-end to use.  See Details.
#' @param nSim Integer; Monte Carlo replicates when `method = "mc"` or the
#'   lookup falls back to MC.
#' @param \dots Currently unused; reserved for future tuning parameters.
#'
#' @returns A numeric vector of normalized expected distances, one entry per
#'   element of `n`.
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
#' @template MRS
#' @references
#' \insertAllCited{}
#' @encoding UTF-8
#' @family tree distances
#' @seealso [`randomTreeDistances`]
#' @export
ExpectedClusteringInfoDistance <- function(n,
                                           method = c("lookup", "enumerate",
                                                      "mc"),
                                           nSim = 1e4L, ...) {
  .ExpectedRandomDistance(n, method = method, nSim = nSim,
                          dist_fun = ClusteringInfoDistance,
                          table_col = "cid_mean")
}

#' @rdname ExpectedClusteringInfoDistance
#' @export
ExpectedPhylogeneticInfoDistance <- function(n,
                                             method = c("lookup", "enumerate",
                                                        "mc"),
                                             nSim = 1e4L, ...) {
  .ExpectedRandomDistance(n, method = method, nSim = nSim,
                          dist_fun = PhylogeneticInfoDistance,
                          table_col = "pid_mean")
}

# Internal: shared back-end for the Expected... functions.
.ExpectedRandomDistance <- function(n, method, nSim, dist_fun, table_col) {
  method <- match.arg(method, c("lookup", "enumerate", "mc"))
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

# Shared chance-adjustment engine used by ACID and APhID.
#
# Computes the normalized distance, then maps it to
# `1 - (D - E) / (1 - E)`.  `E` is either supplied via `expected` (scalar or
# matrix matching the shape of `D`) or looked up via `expected_fun` on a
# per-pair basis.  Mirrors the return shape (numeric / dist / matrix) of
# `dist_fun()`.
.AdjustInfoDistance <- function(tree1, tree2, expected, dist_fun,
                                expected_fun, ...) {
  d <- dist_fun(tree1, tree2, normalize = TRUE)
  d_attr <- attributes(d)

  e <- if (is.null(expected)) {
    .ExpectedFor(tree1, tree2, expected_fun, ...)
  } else {
    expected
  }

  # Chance-adjusted similarity:
  #   identity (D = 0)      -> 1
  #   random pair (D ~= E)  -> ~0
  #   more different than random (D > E) -> negative (a real signal)
  adj <- 1 - d / e
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
