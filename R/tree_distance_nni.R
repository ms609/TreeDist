#' Approximate the NNI distance
#' 
#' An approximation of the NNI distance based on Li et al. 1996.
#' 
#' @param tree1,tree2 Individual trees of class `phylo`
#' 
#' @examples
#'   library('TreeTools')
#'   tree1 <- BalancedTree(7)
#'   tree2 <- PectinateTree(7)
#'   
#'   NNIDist(tree1, tree2)
#' 
#' @references 
#'   \insertRef{Li1996}{TreeDist}
#'   
#' @template MRS
#' @importFrom ape Nnode.phylo
#' @importFrom TreeTools PostorderEdges NTip RenumberTips
#' @export
NNIDist <- function (tree1, tree2) {
  edge <- tree1$edge
  edge1 <- do.call(cbind, 
                   PostorderEdges(edge[, 1], edge[, 2], dim(edge)[1], 
                                  Nnode.phylo(tree1)))
  tree2 <- RenumberTips(tree2, tree1$tip.label)
  edge <- tree2$edge
  edge2 <- do.call(cbind, 
                   PostorderEdges(edge[, 1], edge[, 2], dim(edge)[1], 
                                  Nnode.phylo(tree2)))
  cpp_nni_distance(edge1, edge2, NTip(tree1))
}