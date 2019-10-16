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
#' @importFrom TreeSearch LogTreesMatchingSplit LnUnrooted.int
#' @export
RobinsonFouldsInfoSplits <- function (splits1, splits2, reportMatching = FALSE) {
  GeneralizedRF(splits1, splits2,
                function(splits1, splits2, nSplits1, nSplits2) {
    nTip <- dim(splits1)[1]
    lnUnrooted <- LnUnrooted.int(nTip)
    
    ret <- matrix(0, nSplits1, nSplits2)
    for (i in seq_len(nSplits1)) {
      A1 <- splits1[, i]
      for (j in seq_len(nSplits2)) {
        A2 <- splits2[, j]
      
        if (all(A1 == A2) || all(A1 != A2)) {
          nInSplit <- sum(A1)
          ret[i, j] <- LogTreesMatchingSplit(nInSplit, nTip - nInSplit)
        } else {
          ret[i, j] <- lnUnrooted
        }
      }
    }
    
    # Return:
    -(ret - lnUnrooted) / log(2)
  }, maximize = TRUE, reportMatching)
}

#' @describeIn RobinsonFouldsInfo An inefficient implementation of the 
#' Robinson-Foulds distance, included for use with [`VisualizeMatching`].
#' To generate the RF distance efficiently, use the function 
#' \code{\link{ape}{treedist}}.
#' @export
RobinsonFoulds <- function (tree1, tree2, similarity = FALSE, normalize = FALSE,
                                reportMatching = FALSE) {
  unnormalized <- CalculateTreeDistance(RobinsonFouldsSplits, tree1, tree2, 
                                        reportMatching)
  
  NumberOfSplits <- function (tr) {
    if (class(tr) == 'phylo') {
      tr$Nnode - 2L
    } else {
      vapply(tr, function (thisTree) thisTree$Nnode - 2L, 1L)
    }
  }
  
  if (!similarity) unnormalized <- 
    outer(NumberOfSplits(tree1), NumberOfSplits(tree2), '+')[, , drop = TRUE] -
    unnormalized
  
  # Return:
  NormalizeInfo(unnormalized, tree1, tree2, how = normalize,
                InfoInTree = NumberOfSplits, Combine = `+`)
}

#' @describeIn RobinsonFouldsInfo Calculate Robinson-Foulds distance from splits
#' instead of trees.
#' @inheritParams MutualPhylogeneticInfoSplits
#' @export
RobinsonFouldsSplits <- function (splits1, splits2, reportMatching = FALSE) {
  GeneralizedRF(splits1, splits2,
                function(splits1, splits2, nSplits1, nSplits2) {
                  ret <- matrix(0L, nSplits1, nSplits2)
                  for (i in seq_len(nSplits1)) {
                    A1 <- splits1[, i]
                    for (j in seq_len(nSplits2)) {
                      A2 <- splits2[, j]
                      ret[i, j] <- all(A1 == A2) || all(A1 != A2)
                    }
                  }
                  
                  # Return:
                  ret
                }, maximize = TRUE, reportMatching) * 2L
}
