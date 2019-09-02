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
RobinsonFouldsSplits <- function (splits1, splits2,
                                         reportMatching = FALSE) {
  splitsInCommon <- 
  GeneralizedRF(splits1, splits2,
                function(splits1, splits2, nSplits1, nSplits2) {
    matrix((mapply(function(i, j) {
      A1 <- splits1[, i]
      A2 <- splits2[, j]
      
      all(A1 == A2) || all(A1 != A2)
    },  seq_len(nSplits1), rep(seq_len(nSplits2), each=nSplits1)
    )), nSplits1, nSplits2)
                  
  }, maximize = TRUE, reportMatching)
  dim(splits1)[2] - splitsInCommon + dim(splits2)[2] - splitsInCommon
}
