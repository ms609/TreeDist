#' Robinson-Foulds distances, with adjustments for phylogenetic information
#' content
#' 
#' Calculate the Robinson-Foulds distance, or the equivalent similarity 
#' measure, with options to (i) annotate matched splits; (ii) weight splits
#TODO Update reference to Smith
#' according to their phylogenetic information content (Smith 2020).
#' 
#' Note that if `reportMatching = TRUE`, the `pairScores` attribute returns
#' a logical matrix specifying whether each pair of splits is identical.
#' 
#' `InfoRobinsonFoulds()` calculates the tree similarity or distance by summing 
#' the phylogenetic information content of all splits that are (or are not)
#' identical in both trees.  Consequently, splits that are more likely
#' to be identical by chance alone make a smaller contribution to overall
#' tree distance, because their similarity is less remarkable.
#' 
#' @inheritParams TreeDistance
#' @param similarity Logical specifying whether to report the result as a tree
#' similarity, rather than a difference.
#' 
#' @templateVar returns `RobinsonFoulds()` and `InfoRobinsonFoulds()` return
#' @template distReturn
#' 
#' @section Normalization:
#' 
#' - `RobinsonFoulds()` is normalized against the total number of splits that
#'  are present.
#'  
#' - `InfoRobinsonFoulds()` is normalized against the sum of the phylogenetic 
#' information of all splits in both trees, treated independently.
#'  
#' @references 
#' 
#' \insertRef{Robinson1981}{TreeDist}
#' 
#' \insertRef{Steel2006}{TreeDist}
#' 
#' \insertRef{SmithDist}{TreeDist}
#' 
#' @examples
#'  # For BalancedTree, PectinateTree, as.phylo:
#' library('TreeTools', quietly = TRUE, warn.conflicts = FALSE)
#' balanced7 <- BalancedTree(7)
#' pectinate7 <- PectinateTree(7)
#' RobinsonFoulds(balanced7, pectinate7)
#' RobinsonFoulds(balanced7, pectinate7, normalize = TRUE)
#' VisualizeMatching(RobinsonFouldsMatching, balanced7, pectinate7)
#' 
#' InfoRobinsonFoulds(balanced7, pectinate7)
#' VisualizeMatching(InfoRobinsonFoulds, balanced7, pectinate7)
#' @template MRS
#' 
#' 
#' @family tree distances
#' @seealso Display paired splits: [`VisualizeMatching()`]
#' 
#' @export
#' @name Robinson-Foulds

#' @aliases RobinsonFouldsInfo
#' @rdname Robinson-Foulds
InfoRobinsonFoulds <- function (tree1, tree2 = tree1, similarity = FALSE,
                                normalize = FALSE, reportMatching = FALSE) {
  unnormalized <- CalculateTreeDistance(InfoRobinsonFouldsSplits, tree1, tree2, 
                                        reportMatching) * 2
  
  if (!similarity) unnormalized <- 
      outer(SplitwiseInfo(tree1), SplitwiseInfo(tree2), '+')[, , drop = TRUE] -
      unnormalized
  
  # Return:
  NormalizeInfo(unnormalized, tree1, tree2, how = normalize,
                InfoInTree = SplitwiseInfo, Combine = '+')
}

#' @export
RobinsonFouldsInfo <- InfoRobinsonFoulds

#' @rdname Robinson-Foulds
#' @inheritParams SharedPhylogeneticInfoSplits
#' @export
InfoRobinsonFouldsSplits <- function (splits1, splits2, 
                                      nTip = attr(splits1, 'nTip'),
                                      reportMatching = FALSE) {
  
  GeneralizedRF(splits1, splits2, nTip, cpp_robinson_foulds_info,
                maximize = FALSE, reportMatching = reportMatching)
}

#' @rdname Robinson-Foulds
#' @importFrom TreeTools NSplits
#' @export
RobinsonFoulds <- function (tree1, tree2 = tree1, similarity = FALSE,
                            normalize = FALSE, reportMatching = FALSE) {
  unnormalized <- CalculateTreeDistance(RobinsonFouldsSplits, tree1, tree2, 
                                        reportMatching)
  
  if (similarity) unnormalized <- 
    outer(NSplits(tree1), NSplits(tree2), '+')[, , drop = TRUE] -
    unnormalized
  
  # Return:
  NormalizeInfo(unnormalized, tree1, tree2, how = normalize,
                InfoInTree = NSplits, Combine = `+`)
}

#' @describeIn Robinson-Foulds Matched splits, intended for use with 
#' [`VisualizeMatching()`].
#' @param \dots Not used.
#' @importFrom TreeTools NSplits
#' @export
RobinsonFouldsMatching <- function (tree1, tree2 = tree1, similarity = FALSE,
                                    normalize = FALSE, ...) {
  ret <- CalculateTreeDistance(RobinsonFouldsSplits, tree1, tree2,
                                        reportMatching = TRUE)

  ret <- outer(NSplits(tree1), NSplits(tree2), '+')[, , drop = TRUE] -
    ret
  
  attr(ret, 'pairScores') <- !attr(ret, 'pairScores')
  
  # Return:
  ret
}

#' @rdname Robinson-Foulds
#' @inheritParams SharedPhylogeneticInfoSplits
#' @export
RobinsonFouldsSplits <- function (splits1, splits2,
                                  nTip = attr(splits1, 'nTip'),
                                  reportMatching = FALSE) {
  GeneralizedRF(splits1, splits2, nTip, cpp_robinson_foulds_distance,
                maximize = FALSE, reportMatching = reportMatching)
}
