#' Nye _et al_. (2006) tree comparison
#' 
#' Implements the tree comparison metric of Nye _et al_. (2006).
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
#' (tree1, tree2, k = 1, arboreal = FALSE)}.
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
#'  from zero to one by dividing by the maximum value possible for any pair 
#'  of trees with  _n_ leaves, $n - 3$.
#' You may wish to normalize instead against the number of splits present
#' in the smaller tree, which represents the maximum value possible for a pair
#' of trees with the specified topologies (`normalize = pmin.int`), or the
#' number of splits in the most resolved tree (`normalize = pmax.int`). 
#' 
#' If `normalize = TRUE` and `similarity = FALSE`, then results will be rescaled
#' from zero to one by dividing by the total number of splits in the pair
#' of trees being considered.
#' 
#' @examples 
#' library('TreeTools')
#' NyeTreeSimilarity(BalancedTree(8), PectinateTree(8))
#' NyeTreeSimilarity(lapply(0:5, as.phylo, nTip = 8), PectinateTree(8))
#' NyeTreeSimilarity(lapply(0:5, as.phylo, nTip = 8), similarity = FALSE)
#' 
#' @template distReturn
#' 
#' @references \insertRef{Nye2006}{TreeDist}
#' @family tree distances
#' 
#' @template MRS
#' @importFrom TreeTools NSplits TipLabels
#' @export
NyeTreeSimilarity <- function (tree1, tree2 = tree1, similarity = TRUE,
                               normalize = FALSE, normalizeMax = TRUE,
                               reportMatching = FALSE) {
  
  unnormalized <- CalculateTreeDistance(NyeSplitSimilarity, tree1, tree2, 
                                        reportMatching)
  if (similarity) {
    MaxSplits <- function (tree) {
      if (inherits(tree, c('phylo', 'Splits'))) {
        length(TipLabels(tree)) - 3L
      } else if (mode(tree) == 'list') {
        vapply(tree, MaxSplits, integer(1L))
      } else {
        stop ("Unrecognized tree format")
      }
    }
    
    # Return:
    NormalizeInfo(unnormalized, tree1, tree2, how = normalize,
                  InfoInTree = if (normalizeMax) MaxSplits else NSplits,
                  Combine = pmax.int)
  } else {
    unnormalized <- outer(NSplits(tree1), NSplits(tree2), '+')[, , drop=TRUE] -
                    (2 * unnormalized)
    
    # Return:
    NormalizeInfo(unnormalized, tree1, tree2, how = normalize,
                  InfoInTree = NSplits, Combine = '+')
  }
}

#' @rdname NyeTreeSimilarity
#' @inheritParams SharedPhylogeneticInfoSplits
#' @export
NyeSplitSimilarity <- function (splits1, splits2, 
                                nTip = attr(splits1, 'nTip'),
                                reportMatching = FALSE) {
  GeneralizedRF(splits1, splits2, nTip, cpp_jaccard_similarity, k = 1L,
                arboreal = FALSE, maximize = TRUE,
                reportMatching = reportMatching)
}

#' Jaccard-Robinson-Foulds metric
#' 
#' Calculate the 
#' [Jaccard-Robinson-Foulds metric](https://ms609.github.io/TreeDist/articles/Generalized-RF.html#jaccard-robinson-foulds-metric)
#' (B&ouml;cker _et al_. 2013).
#' 
#' In short, the Jaccard-Robinson-Foulds metric is a generalized Robinson-Foulds
#' metric: it finds the optimal matching that pairs each split in one tree with
#' a similar split in the second.
#' Matchings are scored according to the size of the largest split that is 
#' consistent with both of them, normalized against the Jaccard index, and 
#' raised to an arbitrary exponent. 
#' A more detailed explanation is availble in the 
#' [vignettes](https://ms609.github.io/TreeDist/articles/Generalized-RF.html#jaccard-robinson-foulds-metric).
#' 
#' By default, arboreal matchings are enforced. 
#' 
#' Note that the settings `k = 1, arboreal = FALSE, similarity = TRUE`
#' give the similarity metric of Nye _et al_. (2006); a slightly faster
#' implementation of this metric is available as [`NyeTreeSimilarity`].
#' 
#' The examples section below details how to visualize matchings with 
#' non-default parameter values. 
#' 
#' @inheritParams RobinsonFoulds
#' @param k An arbitrary exponent to which to raise the Jaccard index.
#' Integer values greater than one are anticipated by B&ouml;cker _et al_.
#' The Nye _et al_. metric uses `k = 1`.
#' As k &rarr; &infin;, the metric converges to the Robinson-Foulds metric.
#' @param arboreal Logical specifying whether to enforce arboreal matches, by
#' assigning zero similarity to contradictory pairs of splits on an 
#' _ad hoc_ basis.
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
#' JaccardRobinsonFoulds(tree1, tree2, k = 2, arboreal = FALSE)
#' JaccardRobinsonFoulds(tree1, tree2, k = 2, arboreal = TRUE)
#' 
#' JRF2 <- function (tree1, tree2, ...) 
#'   JaccardRobinsonFoulds(tree1, tree2, k = 2, arboreal= TRUE, ...)
#'   
#' VisualizeMatching(JRF2, tree1, tree2, matchZeros = FALSE)
#' 
#' @family tree distances
#' 
#' @template MRS
#' @importFrom TreeTools NSplits
#' @export
JaccardRobinsonFoulds <- function (tree1, tree2 = tree1, k = 1L, 
                                   arboreal = TRUE, similarity = FALSE,
                                   normalize = FALSE, reportMatching = FALSE) {
  unnormalized <- CalculateTreeDistance(JaccardSplitSimilarity, tree1, tree2, 
                                        k = k, arboreal = arboreal, 
                                        reportMatching = reportMatching) * 2L
  if (!similarity) unnormalized <- 
      outer(NSplits(tree1), NSplits(tree2), '+')[, , drop=TRUE] - unnormalized

  NormalizeInfo(unnormalized, tree1, tree2, how = normalize,
                InfoInTree = NSplits, Combine = '+')
}

#' @rdname JaccardRobinsonFoulds
#' @inheritParams SharedPhylogeneticInfoSplits
#' @export
JaccardSplitSimilarity <- function (splits1, splits2,
                                    nTip = attr(splits1, 'nTip'),
                                    k = 1L, arboreal = TRUE,
                                    reportMatching = FALSE) {
  GeneralizedRF(splits1, splits2, nTip, cpp_jaccard_similarity, k = k,
                arboreal = arboreal, maximize = TRUE,
                reportMatching = reportMatching)
}
