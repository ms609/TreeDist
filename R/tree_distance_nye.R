#' Nye _et al_. (2006) tree comparison
#' 
#' Implements the tree comparison metric of Nye _et al_. (2006).
#' In short, this finds the optimal matching that pairs each branch from
#' one tree with a branch in the second, where matchings are scored according to
#' the size of the largest bipartition that is consistent with both of them,
#' normalized against the Jaccard index.
#' 
#' @inheritParams MutualPhylogeneticInfo
#' 
#' @section Normalization:
#' 
#' If `normalize = TRUE`, then results will be rescaled from zero to one by
#' dividing by the maximum possible value for trees of the given topologies,
#' which is four less than the total number of nodes in both trees. 
#' You may wish to normalize instead against the maximum number of nodes present
#' in a pair of trees with _n_ terminals, by specifying `normalize = <some number>`.
#' 
#' @references \insertRef{Nye2006}{TreeDist}
#' @family tree distances
#' 
#' @author Martin R. Smith
#' @export
NyeTreeSimilarity <- function (tree1, tree2, normalize = FALSE,
                               reportMatching = FALSE) {
  unnormalized <- CalculateTreeDistance(NyeSplitSimilarity, tree1, tree2, 
                                        reportMatching)
  NyeInfoCounter <- function (tr) {
    if (class(tr) == 'phylo') {
      tr$Nnode - 2L
    } else {
      vapply(tr, function (thisTree) thisTree$Nnode - 2L, 1L)
    }
  }
  
  NormalizeInfo(unnormalized, tree1, tree2, how = normalize,
                InfoInTree = NyeInfoCounter, Combine = pmax)
}

#' @describeIn NyeTreeSimilarity Calculate tree similarity from splits 
#' instead of trees.
#' @inheritParams MutualPhylogeneticInfoSplits
#' @export
NyeSplitSimilarity <- function (splits1, splits2, normalize = TRUE,
                                reportMatching = FALSE) {
  
  dimSplits1 <- dim(splits1)
  dimSplits2 <- dim(splits2)
  nSplits1 <- dimSplits1[2]
  nSplits2 <- dimSplits2[2]
  if (nSplits1 == 0 || nSplits2 == 0) return (0L)
  nTerminals <- dimSplits1[1]
  if (dimSplits2[1] != nTerminals) {
    stop("Split rows must bear identical labels")
  }
  
  taxonNames1 <- rownames(splits1)
  taxonNames2 <- rownames(splits2)
  
  if (!is.null(taxonNames2)) {
    splits2 <- unname(splits2[taxonNames1, , drop=FALSE])
    splits1 <- unname(splits1) # split2[split1] faster without names
  }
  
  Ars <- function (pir, pjs) {
    sum(pir[pjs]) / sum(pir | pjs)
  }
  
  pairScores <- matrix((mapply(function(i, j) {
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
  
  # Return:
  TreeDistanceReturn(pairScores, maximize = TRUE, reportMatching, 
                     splits1, splits2, taxonNames1)
}
