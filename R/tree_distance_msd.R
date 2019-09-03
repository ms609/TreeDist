#' Matching Split Distance
#' 
#' Implements the Matching Split Distance for unrooted binary phylogenetic 
#' trees of Bogdanowicz and Giaro (2012).
#' 
#' @inheritParams MutualPhylogeneticInfo
#' 
#' @section Normalization:
#' 
#' A normalization value or function must be provided in order to return a
#' normalized value.  If you are aware of a generalised formula, please
#' let me know by
#' \href{https://github.com/ms609/TreeDist/issues/new}{creating a GitHub issue}
#' so that it can be implemented.
#'  
#' @references \insertRef{Bogdanowicz2012}{TreeDist}
#' @family tree distances
#' 
#' @author Martin R. Smith
#' @importFrom TreeSearch LnUnrooted.int
#' @export
MatchingSplitDistance <- function (tree1, tree2, normalize = FALSE,
                                   reportMatching = FALSE) {
  unnormalized <- CalculateTreeDistance(MatchingSplitDistanceSplits, tree1, tree2, 
                                        reportMatching)
  
  # Return:
  NormalizeInfo(unnormalized, tree1, tree2, how = normalize,
                InfoInTree = function (X) stop("Please specify a function to generate a normalizing constant"),
                Combine = max)
}

#' @describeIn MatchingSplitDistance Calculate Matching Split Distance from splits instead of trees.
#' @inheritParams MutualPhylogeneticInfoSplits
#' @export
MatchingSplitDistanceSplits <- function (splits1, splits2, normalize = TRUE, 
                                         reportMatching = FALSE) {
  GeneralizedRF(splits1, splits2, 
                function(splits1, splits2, nSplits1, nSplits2) {
    SymmetricDifference <- function (A, B) {
      (A & !B) | (!A & B)
    }
    
    matrix((mapply(function(i, j) {
      A1 <- splits1[, i]
      A2 <- splits2[, j]
      B1 <- !A1
      B2 <- !A2
      
      # Long-winded way:
      # min(
      #   sum(SymmetricDifference(A1, A2), SymmetricDifference(B1, B2)),
      #   sum(SymmetricDifference(A1, B2), SymmetricDifference(B1, A2))
      # ) / 2L
      # But SD(A1, A2) == SD(B1, B2) and SD(A1, B2) == SD(B1, A2), so:
      
      min(sum(SymmetricDifference(A1, A2)), sum(SymmetricDifference(A1, B2)))
    },  seq_len(nSplits1), rep(seq_len(nSplits2), each=nSplits1)
    )), nSplits1, nSplits2)
  }, maximize = FALSE, reportMatching)
}
