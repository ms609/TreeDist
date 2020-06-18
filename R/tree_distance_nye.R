#' Nye _et al_. (2006) tree comparison
#' 
#' Implements the 
#' [Generalized Robinson-Foulds](https://ms609.github.io/TreeDist/articles/Robinson-Foulds.html#generalized-robinson-foulds-distances).
#' tree comparison metric of Nye _et al_. (2006).
#' In short, this finds the optimal matching that pairs each branch from
#' one tree with a branch in the second, where matchings are scored according to
#' the size of the largest split that is consistent with both of them,
#' normalized against the Jaccard index.
#' A more detailed account is available in the 
#' [vignettes](https://ms609.github.io/TreeDist/articles/Generalized-RF.html#nye-et-al--tree-similarity-metric).
#' 
#' The measure is defined as a similarity score.  If `similarity = FALSE`, the
#' similarity score will be converted into a distance by doubling it and
#' subtracting it from the number of splits present in both trees.
#' This ensures consistency with `JaccardRobinsonFoulds`.
#' 
#' Note that `NyeTreeSimilarity(tree1, tree2)` is equivalent to, but 
#' slightly faster than, \code{\link{JaccardRobinsonFoulds}
#' (tree1, tree2, k = 1, allowConflict = TRUE)}.
#'  
#' @inheritParams RobinsonFoulds
#' @param normalizeMax When calculating similarity, normalize against the 
#' maximum number of splits that could have been present (`TRUE`),
#'  or the number of splits that were actually observed (`FALSE`)?  
#' Defaults to the number of splits in the better-resolved tree; set
#'  `normalize = pmin.int` to use the number of splits in the less resolved
#'  tree.
#' 
#' @section Normalization:
#' 
#' If `normalize = TRUE` and `similarity = TRUE`, then results will be rescaled
#' from zero to one by dividing by the mean number of splits in the two trees
#' being compared.
#'  
#' You may wish to normalize instead against the number of splits present
#' in the smaller tree, which represents the maximum value possible for a pair
#' of trees with the specified topologies (`normalize = pmin.int`); the
#' number of splits in the most resolved tree (`normalize = pmax.int`);
#' or the maximum value possible for any pair of trees with  _n_ leaves, 
#' _n_ - 3 (`normalize = TreeTools::NTip(tree1) - 3L`).
#' 
#' If `normalize = TRUE` and `similarity = FALSE`, then results will be rescaled
#' from zero to one by dividing by the total number of splits in the pair
#' of trees being considered.
#' 
#' @examples 
#' library('TreeTools')
#' NyeTreeSimilarity(BalancedTree(8), PectinateTree(8))
#' VisualizeMatching(NyeTreeSimilarity ,BalancedTree(8), PectinateTree(8))

#' NyeTreeSimilarity(as.phylo(0:5, nTip = 8), PectinateTree(8))
#' NyeTreeSimilarity(as.phylo(0:5, nTip = 8), similarity = FALSE)
#' 
#' @template distReturn
#' 
#' @references \insertRef{Nye2006}{TreeDist}
#' @family tree distances
#' 
#' @template MRS
#' @importFrom TreeTools NSplits
#' @export
NyeTreeSimilarity <- function (tree1, tree2 = tree1, similarity = TRUE,
                               normalize = FALSE,
                               normalizeMax = !is.logical(normalize),
                               reportMatching = FALSE) {
  
  unnormalized <- CalculateTreeDistance(NyeSplitSimilarity, tree1, tree2, 
                                        reportMatching)
  if (similarity) {
    # Return:
    NormalizeInfo(unnormalized, tree1, tree2, how = normalize,
                  InfoInTree = if (normalizeMax) SplitsInBinaryTree else NSplits,
                  Combine = .MeanOfTwo)
  } else {
    unnormalized <- outer(NSplits(tree1), NSplits(tree2), '+')[, , drop = TRUE] -
                    (unnormalized + unnormalized)
    
    # Return:
    NormalizeInfo(unnormalized, tree1, tree2, how = normalize,
                  InfoInTree = NSplits, Combine = '+')
  }
}

.MeanOfTwo <- function (x, y) (x + y) / 2L

#' @rdname NyeTreeSimilarity
#' @inheritParams SharedPhylogeneticInfoSplits
#' @export
NyeSplitSimilarity <- function (splits1, splits2, 
                                nTip = attr(splits1, 'nTip'),
                                reportMatching = FALSE) {
  GeneralizedRF(splits1, splits2, nTip, cpp_jaccard_similarity, k = 1L,
                allowConflict = TRUE, maximize = TRUE,
                reportMatching = reportMatching)
}

#' Jaccard-Robinson-Foulds metric
#' 
#' Calculate the 
#' [Jaccard-Robinson-Foulds metric](https://ms609.github.io/TreeDist/articles/Generalized-RF.html#jaccard-robinson-foulds-metric)
#' (B&ouml;cker _et al_. 2013), a 
#' [Generalized Robinson-Foulds metric](https://ms609.github.io/TreeDist/articles/Robinson-Foulds.html#generalized-robinson-foulds-distances).
#' 
#' In short, the Jaccard-Robinson-Foulds metric is a generalized Robinson-Foulds
#' metric: it finds the optimal matching that pairs each split in one tree with
#' a similar split in the second.
#' Matchings are scored according to the size of the largest split that is 
#' consistent with both of them, normalized against the Jaccard index, and 
#' raised to an arbitrary exponent. 
#' A more detailed explanation is provided in the 
#' [vignettes](https://ms609.github.io/TreeDist/articles/Generalized-RF.html#jaccard-robinson-foulds-metric).
#' 
#' By default, conflicting splits may be paired. 
#' 
#' Note that the settings `k = 1, allowConflict = TRUE, similarity = TRUE`
#' give the similarity metric of Nye _et al_. (2006); a slightly faster
#' implementation of this metric is available as [`NyeTreeSimilarity()`].
#' 
#' The examples section below details how to visualize matchings with 
#' non-default parameter values.
#' 
#' @inheritParams RobinsonFoulds
#' @param k An arbitrary exponent to which to raise the Jaccard index.
#' Integer values greater than one are anticipated by B&ouml;cker _et al_.
#' The Nye _et al_. metric uses `k = 1`.
#' As k increses towards infinity, the metric converges to the Robinson-Foulds
#' metric.
#' @param allowConflict Logical specifying whether to allow conflicting splits
#' to be paired. If `FALSE`, such pairings will be allocated a similarity
#' score of zero.
#' 
#' @section Normalization:
#' 
#' If `normalize = TRUE`, then results will be rescaled from zero to one by
#' dividing by the maximum possible value for trees of the given topologies,
#' which is equal to the sum of the number of splits in each tree. 
#' You may wish to normalize instead against the maximum number of splits
#' present in a pair of trees with _n_ leaves, by specifying 
#' `normalize = n - 3`.
#' 
#' @template distReturn
#' 
#' @references 
#' 
#' - \insertRef{Nye2006}{TreeDist}
#' 
#' - \insertRef{Bocker2013}{TreeDist}
#' 
#' @examples 
#' set.seed(2)
#' tree1 <- ape::rtree(10)
#' tree2 <- ape::rtree(10)
#' JaccardRobinsonFoulds(tree1, tree2, k = 2, allowConflict = FALSE)
#' JaccardRobinsonFoulds(tree1, tree2, k = 2, allowConflict = TRUE)
#' 
#' JRF2 <- function (tree1, tree2, ...) 
#'   JaccardRobinsonFoulds(tree1, tree2, k = 2, allowConflict = FALSE, ...)
#'   
#' VisualizeMatching(JRF2, tree1, tree2, matchZeros = FALSE)
#' @template MRS
#' 
#' @family tree distances
#' 
#' @encoding UTF-8
#' @importFrom TreeTools NSplits
#' @export
JaccardRobinsonFoulds <- function (tree1, tree2 = tree1, k = 1L, 
                                   allowConflict = TRUE, similarity = FALSE,
                                   normalize = FALSE, reportMatching = FALSE) {
  unnormalized <- CalculateTreeDistance(JaccardSplitSimilarity, tree1, tree2, 
                                        k = k, allowConflict = allowConflict, 
                                        reportMatching = reportMatching) * 2L
  if (!similarity) unnormalized <- 
      outer(NSplits(tree1), NSplits(tree2), '+')[, , drop = TRUE] - unnormalized
#TODO make normalization match Nye Et Al.
  NormalizeInfo(unnormalized, tree1, tree2, how = normalize,
                InfoInTree = NSplits, Combine = '+')
}

#' @rdname JaccardRobinsonFoulds
#' @inheritParams SharedPhylogeneticInfoSplits
#' @export
JaccardSplitSimilarity <- function (splits1, splits2,
                                    nTip = attr(splits1, 'nTip'),
                                    k = 1L, allowConflict = TRUE,
                                    reportMatching = FALSE) {
  GeneralizedRF(splits1, splits2, nTip, cpp_jaccard_similarity, k = k,
                allowConflict = allowConflict, maximize = TRUE,
                reportMatching = reportMatching)
}
