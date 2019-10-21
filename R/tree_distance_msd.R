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
#' @importFrom TreeTools LnUnrooted.int
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
#' @useDynLib TreeDist, .registration = TRUE
#' @export
MatchingSplitDistanceSplits <- function (splits1, splits2, normalize = TRUE, 
                                         reportMatching = FALSE) {
  CGRF(splits1, splits2, cpp_matching_split_distance, maximize = FALSE,
       reportMatching = reportMatching)
}
