#' Median of a set of trees
#' 
#' Provides a single binary tree that represents the geometric median -- 
#' an 'average' -- of a pool of tree topologies.
#' 
#' The geometric median is the tree that exhibits the shortest average distance
#' from each other tree topology in the set.
#' It represents an 'average' of a set of trees, though note that an unsampled 
#' tree may be closer to the geometric 'centre of gravity' of the input set --
#' such a tree would not be considered.
#' 
#' The result will depend on the metric chosen to calculate distances between
#' tree topologies. In the absence of a natural metric of tree topologies, 
#' the default choice is [`ClusteringInfoDistance()`] -- which discards
#' branch length information.
#' If specifying a different function, be sure that it returns a difference,
#' rather than a similarity.
#' 
#' 
#' @param x Object of class `multiPhylo` containing phylogenetic trees.
#' @param na.rm,\dots Unused; included for consistency with default function..
#' @param Distance Function to calculate distances between each pair
#' of trees in `x`.
#' @param index Logical: if `TRUE`, return the index of the median tree(s);
#' if `FALSE`, return the tree itself.
#' @param breakTies Logical: if `TRUE`, return a single tree with the minimum
#'  score; if `FALSE`, return all tied trees.
#' 
#' @return `median.multiPhylo()` returns an object of class `phylo` 
#' corresponding to the geometric median of a set of trees:
#' that is, the tree whose average distance from all other trees in the set
#' is lowest.  
#' If multiple trees tie in their average distance, the first will be returned,
#' unless `breakTies = TRUE`, in which case an object of class `multiPhylo`
#' containing all such trees will be returned.
#' 
#' @examples
#' library('TreeTools')
#' tenTrees <- as.phylo(1:10, nTip = 8)
#' 
#' # Default settings:
#' median(tenTrees)
#' 
#' # The Robinson-Foulds distance will create ties:
#' median(tenTrees, Distance = RobinsonFoulds, breakTies = FALSE)
#' 
#' # Be sure to use a distance function, rather than a similarity:
#' NyeDistance <- function (...) NyeSimilarity(..., similarity = FALSE)
#' median(tenTrees, Distance = NyeDistance)
#' 
#' # To analyse a list of trees that is not of class multiPhylo:
#' treeList <- lapply(1:10, as.phylo, nTip = 8)
#' median(structure(treeList, class = 'multiPhylo'))
#' 
#' @seealso Consensus methods:
#'   [`ape::consensus`], 
#'   [`TreeTools::ConsensusWithout`]
#'   
#' @template MRS
#' @importFrom stats median
#' @export
median.multiPhylo <- function (x, na.rm = FALSE, 
                               Distance = ClusteringInfoDistance,
                               index = FALSE,
                               breakTies = TRUE, ...) {
  distances <- colSums(Distance(x))
  
  # Return:
  if (breakTies) {
    if (index) which.min(distances) else x[[which.min(distances)]]
  } else {
    chosen <- distances == min(distances)
    if (index) {
      which(chosen) 
    } else {
      if (sum(chosen) == 1L) x[[which(chosen)]] else x[chosen]
    }
  }
}
