#' Adjusted Clustering Information Distance
#'
#' Chance-adjusted variant of [`ClusteringInfoDistance()`]
#' \insertCite{SmithDist}{TreeDist} that rescales the normalized distance
#' relative to the expected distance between a random pair of `n`-leaf binary
#' trees, so that 1 denotes identical trees and 0 denotes "no more
#' information in common than a random pair".
#'
#' For the normalized clustering information distance `D_norm` (in `[0, 1]`)
#' with expected value `E[D_norm | n]` over random `n`-leaf binary tree pairs,
#' the adjusted score is
#'
#'   `D_adj = 1 - D_norm / E[D_norm | n]`,
#'
#' so that identical trees score 1, a pair drawn from the random distribution
#' scores ~0 on average, and pairs that are more different than expected by
#' chance score below 0.
#'
#' `E[D_norm | n]` is read from the shipped reference table
#' [`randomTreeDistances`] via [`ExpectedClusteringInfoDistance()`]; for tip
#' counts beyond the shipped range an on-the-fly Monte Carlo fallback is used
#' (and a message emitted).  When `tree1` and `tree2` mix trees with different
#' tip counts, the effective `n` for each pair is the number of tips shared
#' between the two trees (since [`ClusteringInfoDistance()`] drops
#' unshared tips before comparison).
#'
#' @template tree12ListParams
#' @param expected Optional numeric scalar or matrix providing
#'   `E[D_norm | n]` directly; if `NULL` (the default) it is looked up via
#'   [`ExpectedClusteringInfoDistance()`].  Useful for benchmarking against
#'   alternative null models.
#' @param null Character; null model used to compute the expected distance.
#'   `"size"` (default) conditions only on the number of tips.  `"shape"`
#'   additionally fixes the unlabelled tree shape, randomising only the tip
#'   labels; see [`ExpectedClusteringInfoDistance()`] for details.  When
#'   `null = "shape"`, both `tree1` and `tree2` must be single `phylo`
#'   objects.
#' @param \dots Additional arguments passed to
#'   [`ExpectedClusteringInfoDistance()`] (e.g. `method`, `nSim`).
#'
#' @returns A numeric vector, `dist` object, or matrix of adjusted clustering
#'   information similarities, matching the structure returned by
#'   [`ClusteringInfoDistance()`].  Values near 1 indicate identical trees;
#'   values near 0 indicate pairs no more similar than two random binary
#'   trees with the same number of tips; negative values indicate pairs more
#'   different than chance.  When `E[D_norm | n] == 0` (a degenerate case
#'   that should not arise in practice), the result is `NA` with a warning.
#'
#' @examples
#' library("TreeTools", quietly = TRUE)
#' t1 <- BalancedTree(8)
#' t2 <- PectinateTree(8)
#' AdjustedClusteringInfoDistance(t1, t1) # 1
#' AdjustedClusteringInfoDistance(t1, t2)
#'
#' # Shape-conditioned null (requires shipped lookup or shape-mc)
#' \dontrun{
#' AdjustedClusteringInfoDistance(t1, t2, null = "shape")
#' }
#'
#' @template MRS
#' @references
#' \insertAllCited{}
#' @encoding UTF-8
#' @family tree distances
#' @seealso [`ClusteringInfoDistance()`], [`ExpectedClusteringInfoDistance()`]
#' @export
AdjustedClusteringInfoDistance <- function(tree1, tree2 = NULL,
                                           expected = NULL,
                                           null = c("size", "shape"),
                                           ...) {
  null <- match.arg(null)
  .AdjustInfoDistance(tree1, tree2, expected = expected,
                      dist_fun = ClusteringInfoDistance,
                      expected_fun = ExpectedClusteringInfoDistance,
                      null = null, ...)
}
