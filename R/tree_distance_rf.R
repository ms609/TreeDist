#' (Information-adjusted) Robinson-Foulds distance
#' 
#' Calculates the Robinson-Foulds distance, or the equivalent similarity 
#' measure, optionally annotating matched partitions and weighting partitions
#' according to their phylogenetic information content.
#' 
#' `RobinsonFoulds` is an inefficient implementation of the Robinson-Foulds 
#' distance, included for use with [`VisualizeMatching`]. 
#' To generate the RF distance efficiently,
#' use the function \code{\link{ape}{treedist}}.
#' 
#' Note that if `reportMatching = TRUE`, the `pairScores` attribute returns
#' a logical matrix specifying whether each pair of partitions is identical.
#' 
#' `RobinsonFouldsInfo` calculates the tree similarity or distance by summing 
#' the phylogenetic information content of all partitions that are (or are not)
#' identical in both trees.  Consequently, partitions that are more likely
#' to be identical by chance alone make a smaller contribution to overall
#' tree distance, because their similarity is less remarkable.
#' 
#' @inheritParams MutualPhylogeneticInfo
#' @param similarity Logical specifying whether to report the result as a tree
#' similarity, rather than a difference.
#' 
#' @section Normalization:
#' 
#' - `RobinsonFoulds` is normalized against the total number of partitions that
#'  are present.
#'  
#' - `RobinsonFouldsInfo` is normalized against the sum of the phylogenetic 
#' information of all splits in both trees, treated independently.
#'  
#' @references \insertRef{Robinson1981}{TreeDist}
#' @family tree distances
#' 
#' @author Martin R. Smith
#' @export
RobinsonFouldsInfo <- function (tree1, tree2, similarity = FALSE,
                                normalize = FALSE, reportMatching = FALSE) {
  unnormalized <- CalculateTreeDistance(RobinsonFouldsInfoSplits, tree1, tree2, 
                                        reportMatching) * 2
  
  if (!similarity) unnormalized <- 
      outer(PartitionInfo(tree1), PartitionInfo(tree2), '+')[, , drop=TRUE] -
      unnormalized
  
  # Return:
  NormalizeInfo(unnormalized, tree1, tree2, how = normalize,
                InfoInTree = PartitionInfo, Combine = '+')
}

#' @describeIn RobinsonFouldsInfo Calculate information-adjusted Robinson-Foulds
#' distance from splits instead of trees.
#' @inheritParams MutualPhylogeneticInfoSplits
#' @importFrom TreeTools LogTreesMatchingSplit LnUnrooted.int
#' @export
RobinsonFouldsInfoSplits <- function (splits1, splits2, 
                                      nTip = attr(splits1, 'nTip'),
                                      reportMatching = FALSE) {
  
  CGRF(splits1, splits2, nTip, cpp_robinson_foulds_info,
       maximize = FALSE, reportMatching = reportMatching)
}

#' @describeIn RobinsonFouldsInfo Robinson-Foulds distance, with option to
#' report matched splits.
#' @importFrom TreeTools NSplits
#' @export
RobinsonFoulds <- function (tree1, tree2, similarity = FALSE, normalize = FALSE,
                                reportMatching = FALSE) {
  unnormalized <- CalculateTreeDistance(RobinsonFouldsSplits, tree1, tree2, 
                                        reportMatching)
  
  if (similarity) unnormalized <- 
    outer(NSplits(tree1), NSplits(tree2), '+')[, , drop = TRUE] -
    unnormalized
  
  # Return:
  NormalizeInfo(unnormalized, tree1, tree2, how = normalize,
                InfoInTree = NSplits, Combine = `+`)
}

#' @describeIn RobinsonFouldsInfo Matched splits, intended for use with 
#' [`VisualizeMatching`].
#' @importFrom TreeTools NSplits
RobinsonFouldsMatching <- function (tree1, tree2, similarity = FALSE,
                                    normalize = FALSE, ...) {
  ret <- CalculateTreeDistance(RobinsonFouldsSplits, tree1, tree2, 
                                        reportMatching = TRUE)

  ret <- outer(NSplits(tree1), NSplits(tree2), '+')[, , drop = TRUE] -
    ret
  
  attr(ret, 'pairScores') <- !attr(ret, 'pairScores')
  
  # Return:
  ret
}

#' @describeIn RobinsonFouldsInfo Calculate Robinson-Foulds distance from splits
#' instead of trees.
#' @inheritParams MutualPhylogeneticInfoSplits
#' @export
RobinsonFouldsSplits <- function (splits1, splits2,
                                  nTip = attr(splits1, 'nTip'),
                                  reportMatching = FALSE) {
  CGRF(splits1, splits2, nTip, cpp_robinson_foulds_distance,
       maximize = FALSE, reportMatching = reportMatching)
}
