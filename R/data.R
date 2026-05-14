#' Shape-conditioned expected information-theoretic distances
#'
#' Precomputed lookup table used by [`ExpectedClusteringInfoDistance()`] and
#' [`ExpectedPhylogeneticInfoDistance()`] when `null = "shape"` and
#' `method = "shape-lookup"`.  Loaded lazily on first use via
#' [`.ShapeExpectedDistances()`] from `inst/extdata/shapeExpectedDistances.rds`.
#'
#' Each row records the mean and standard deviation of the (normalized)
#' [`ClusteringInfoDistance()`] and [`PhylogeneticInfoDistance()`] for a pair
#' of `n`-leaf trees with unlabelled shapes `shape1` and `shape2` under
#' a random label-shuffle null (tip labels drawn independently and uniformly
#' at random for each shape).  Values are Monte Carlo estimates (1 000
#' replicates per pair) computed by `data-raw/build_shape_lookup.R`.
#' Shape integers are **0-based** (range `0:(NRootedShapes(n) - 1)`).
#'
#' The table covers `n = 4:10`.  For `n > 10`, use `method = "shape-mc"`.
#'
#' @format A named list, one element per `n` (names `"4"`, ..., `"10"`).
#' Each element is a data frame with one row per shape pair (`shape1 <= shape2`)
#' and columns
#' \describe{
#'   \item{`shape1`, `shape2`}{Integer (0-based) shape indices.}
#'   \item{`mean_cid`, `sd_cid`}{Mean and SD of the normalized clustering
#'     information distance over `n_pairs` Monte Carlo replicates.}
#'   \item{`mean_pid`, `sd_pid`}{Mean and SD of the normalized phylogenetic
#'     information distance.}
#'   \item{`n_pairs`}{Number of Monte Carlo replicates used.}
#' }
#'
#' @seealso [`ExpectedClusteringInfoDistance()`], [`.ShapeExpectedDistances()`]
#' @encoding UTF-8
#' @keywords datasets
"shapeExpectedDistances"

#' Expected information-theoretic distances between random tree pairs
#'
#' Precomputed reference table used by [`AdjustedClusteringInfoDistance()`] and
#' [`AdjustedPhylogeneticInfoDistance()`] to chance-adjust their scores.
#'
#' Each row records the mean and standard deviation of the (normalized)
#' [`ClusteringInfoDistance()`] and [`PhylogeneticInfoDistance()`] for a random
#' pair of `n`-leaf binary trees.  For `n = 4`...`7` the values are computed
#' by **exact enumeration** over every ordered pair of labelled rooted binary
#' trees (including self-pairs, matching the independent-draw Monte Carlo
#' protocol); for `n >= 8` they are Monte Carlo estimates drawn with
#' [`TreeTools::RandomTree()`], using the protocol of
#' `TreeDistData::randomTreeDistances` \insertCite{SmithDist}{TreeDist}
#' topped up so that `SE / mean <= 0.001` per cell.  See
#' `data-raw/enumerate_small_n.R` and `data-raw/build_random_tree_distances.R`
#' for the build pipeline.
#'
#' @format A data frame with one row per tip count (`n = 4`...`200`) and
#' columns
#' \describe{
#'   \item{`n`}{Integer; number of leaves.}
#'   \item{`cid_mean`,`cid_sd`}{Mean and standard deviation of the normalized
#'     clustering information distance.}
#'   \item{`pid_mean`,`pid_sd`}{Mean and standard deviation of the normalized
#'     phylogenetic information distance.}
#'   \item{`n_repls_cid`,`n_repls_pid`}{Number of replicate pairs contributing
#'     to the CID and PID estimates respectively.  `Inf` flags an exact
#'     enumeration value (n = 4...7).}
#' }
#'
#' @references
#' \insertAllCited{}
#' @encoding UTF-8
#' @keywords datasets
"randomTreeDistances"
