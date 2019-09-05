#' Nye _et al_. (2006) tree comparison
#' 
#' Implements the tree comparison metric of Nye _et al_. (2006).
#' In short, this finds the optimal matching that pairs each branch from
#' one tree with a branch in the second, where matchings are scored according to
#' the size of the largest bipartition that is consistent with both of them,
#' normalized against the Jaccard index.  This is equivalent to (but 
#' slightly faster than) \code{\link{JaccardRobinsonFoulds}
#' (tree1, tree2, k = 1, arboreal = FALSE)}.
#' 
#' @inheritParams RobinsonFoulds
#' 
#' @section Normalization:
#' 
#' If `normalize = TRUE`, then results will be rescaled from zero to one by
#' dividing by the maximum possible value for trees of the given topologies,
#' which is four less than the total number of nodes in both trees. 
#' You may wish to normalize instead against the maximum number of nodes present
#' in a pair of trees with _n_ terminals, by specifying 
#' `normalize = <some number>`.
#' 
#' @references \insertRef{Nye2006}{TreeDist}
#' @family tree distances
#' 
#' @author Martin R. Smith
#' @importFrom TreeSearch NPartitions
#' @export
NyeTreeSimilarity <- function (tree1, tree2, similarity = TRUE,
                               normalize = FALSE, reportMatching = FALSE) {
  unnormalized <- CalculateTreeDistance(NyeSplitSimilarity, tree1, tree2, 
                                        reportMatching)
  if (!similarity) unnormalized <- 
      outer(NPartitions(tree1), NPartitions(tree2), '+')[, , drop=TRUE] -
      unnormalized
  NormalizeInfo(unnormalized, tree1, tree2, how = normalize,
                InfoInTree = NPartitions, Combine = pmax)
}

#' @describeIn NyeTreeSimilarity Calculate tree similarity from splits 
#' instead of trees.
#' @inheritParams MutualPhylogeneticInfoSplits
#' @export
NyeSplitSimilarity <- function (splits1, splits2, reportMatching = FALSE) {
  GeneralizedRF(splits1, splits2, 
                function(splits1, splits2, nSplits1, nSplits2) {
    Ars <- function (pir, pjs) {
      sum(pir[pjs]) / sum(pir | pjs)
    }
    
    # Return:
    matrix((mapply(function(i, j) {
      splitI0 <- splits1[, i]
      splitJ0 <- splits2[, j]
      splitI1 <- !splitI0
      splitJ1 <- !splitJ0
      
      max(
        min(Ars(splitI0, splitJ0), Ars(splitI1, splitJ1)),
        min(Ars(splitI0, splitJ1), Ars(splitI1, splitJ0))
      )
      
    }, seq_len(nSplits1), rep(seq_len(nSplits2), each=nSplits1)
    )), nSplits1, nSplits2)
  }, maximize = TRUE, reportMatching)
}

#' Jaccard-Robinson-Foulds metric
#' 
#' Implements the Jaccard-Robinson-Foulds metric of B&ouml;cker _et al_. (2013).
#' 
#' In short, this finds the optimal matching that pairs each branch from
#' one tree with a branch in the second, where matchings are scored according to
#' the size of the largest bipartition that is consistent with both of them,
#' normalized against the Jaccard index, and raised to an arbitrary exponent.
#' 
#' By default, arboreal matchings are enforced. 
#' 
#' Note that the settings `k = 1, arboreal = FALSE` give the similarity metric
#' of Nye _et al_. (2006); a slightly faster implementation of this metric is
#' available as `[NyeTreeSimilarity]`.
#' 
#' The examples section details how to visualize matchings with non-default
#' parameter values. 
#' 
#' @inheritParams RobinsonFoulds
#' @param k An arbitrary exponent to which to raise the Jaccard index.
#' Integer values greater than one are anticipated by B&ouml;cker _et al_.
#' The Nye _et al_. metric uses `k = 1`.
#' As k → &infin;, the metric converges to the Robinson-Foulds metric.
#' @param arboreal Logical specifying whether to enforce arboreal matches, by
#' assigning zero similarity to contradictory pairs of partitions on an 
#' _ad hoc_ basis.
#' 
#' @section Normalization:
#' 
#' If `normalize = TRUE`, then results will be rescaled from zero to one by
#' dividing by the maximum possible value for trees of the given topologies,
#' which is four less than the total number of nodes in both trees. 
#' You may wish to normalize instead against the maximum number of nodes present
#' in a pair of trees with _n_ terminals, by specifying 
#' `normalize = <some number>`.
#' 
#' @references 
#' 
#' - \insertRef{Nye2006}{TreeDist}
#' 
#' - \insertRef{Bocker2013}{TreeDist}
#' 
#' @examples {
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
#' }
#' 
#' 
#' @family tree distances
#' 
#' @author Martin R. Smith
#' @importFrom TreeSearch NPartitions
#' @export
JaccardRobinsonFoulds <- function (tree1, tree2, k = 1L, arboreal = TRUE,
                                   similarity = FALSE,
                                   normalize = FALSE, reportMatching = FALSE) {
  unnormalized <- CalculateTreeDistance(JaccardSplitSimilarity, tree1, tree2, 
                                        k = k, arboreal = arboreal, 
                                        reportMatching = reportMatching) * 2L
  if (!similarity) unnormalized <- 
      outer(NPartitions(tree1), NPartitions(tree2), '+')[, , drop=TRUE] -
      unnormalized

  NormalizeInfo(unnormalized, tree1, tree2, how = normalize,
                InfoInTree = NPartitions, Combine = '+')
}

#' @describeIn JaccardRobinsonFoulds Calculate tree similarity from splits 
#' instead of trees.
#' @inheritParams MutualPhylogeneticInfoSplits
#' @export
JaccardSplitSimilarity <- function (splits1, splits2,
                                    k = 1L, arboreal = TRUE,
                                    reportMatching = FALSE) {
  GeneralizedRF(splits1, splits2, 
                function(splits1, splits2, nSplits1, nSplits2) {
    Ars <- function (pir, pjs) {
      sum(pir[pjs]) / sum(pir | pjs)
    }
    
    # Return:
    matrix((mapply(function(i, j) {
      splitI0 <- splits1[, i]
      splitJ0 <- splits2[, j]
      splitI1 <- !splitI0
      splitJ1 <- !splitJ0
      
      if (arboreal && !(
        all(splitI0[splitJ0]) ||
        all(splitI0[splitJ1]) ||
        all(splitI1[splitJ0]) ||
        all(splitI1[splitJ1]))) {
        0
      } else {
        max(
          min(Ars(splitI0, splitJ0), Ars(splitI1, splitJ1)),
          min(Ars(splitI0, splitJ1), Ars(splitI1, splitJ0))
        ) ^ k
      }
    }, seq_len(nSplits1), rep(seq_len(nSplits2), each=nSplits1)
    )), nSplits1, nSplits2)
  }, maximize = TRUE, reportMatching)
}
