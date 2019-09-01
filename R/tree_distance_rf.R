#' Robinson-Foulds Distance
#' 
#' An inefficient implementation of the Robinson-Foulds distance, included
#' for use with [`VisualizeMatching`].  To generate the RF distance efficiently,
#' use the function \code{\link{ape}{treedist}}.
#' 
#' @inheritParams MutualPhylogeneticInfo
#' 
#' @section Normalization:
#' 
#' Normalized against the total number of partitions that are present.
#'  
#' @references \insertRef{Robinson1981}{TreeDist}
#' @family tree distances
#' 
#' @author Martin R. Smith
#' @export
RobinsonFoulds <- function (tree1, tree2, normalize = FALSE,
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
  
  # Return:
  NormalizeInfo(unnormalized, tree1, tree2, how = normalize,
                InfoInTree = NumberOfSplits, Combine = `+`)
}

#' @describeIn RobinsonFoulds Calculate Robinson-Foulds distance from splits
#' instead of trees.
#' @inheritParams MutualPhylogeneticInfoSplits
#' @importFrom clue solve_LSAP
#' @export
RobinsonFouldsSplits <- function (splits1, splits2, normalize = TRUE, 
                                         reportMatching = FALSE) {
  
  dimSplits1 <- dim(splits1)
  dimSplits2 <- dim(splits2)
  nTerminals <- dimSplits1[1]
  if (dimSplits2[1] != nTerminals) {
    stop("Split rows must bear identical labels")
  }
  
  swapSplits <- (dimSplits1[2] > dimSplits2[2])
  if (swapSplits) {
    # solve_LDAP expects splits1 to be no larger than splits2
    tmp <- splits1
    splits1 <- splits2
    splits2 <- tmp
    
    tmp <- dimSplits1
    dimSplits1 <- dimSplits2
    dimSplits2 <- tmp
    
    remove(tmp)
  }
  
  taxonNames1 <- rownames(splits1)
  taxonNames2 <- rownames(splits2)
  
  if (!is.null(taxonNames2)) {
    splits2 <- unname(splits2[taxonNames1, , drop=FALSE])
    splits1 <- unname(splits1) # split2[split1] faster without names
  }
  
  nSplits1 <- dimSplits1[2]
  nSplits2 <- dimSplits2[2]
  if (nSplits1 == 0) return (0)
  
  pairScores <- 1 - matrix((mapply(function(i, j) {
    A1 <- splits1[, i]
    A2 <- splits2[, j]
    B1 <- !A1
    B2 <- !A2
    
    all(A1 == A2) || all(A1 == B2)
  },  seq_len(nSplits1), rep(seq_len(nSplits2), each=nSplits1)
  )), nSplits1, nSplits2)
  
  # Return:
  ret <- TreeDistanceReturn(pairScores, reportMatching, swapSplits,
                     taxonNames1)
  
  Score = function (pairScores, matchedSplits)
    sum(dim(pairScores), -matchedSplits, -matchedSplits)
}
