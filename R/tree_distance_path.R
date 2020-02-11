#' Path distance
#' 
#' Calculate the path distance between trees.
#' 
#' This function is a wrapper for the function 
#' \code{\link[phangorn:treedist]{path.dist()}} in the phangorn package.
#' It pre-processes trees to ensure that their internal representation does
#' not cause the `path.dist()` function to crash R.
#' 
#' @template tree12Params
#' 
#' @return `PathDist()` returns a vector or matrix of SPR distances between
#'  trees.
#' 
#' @examples
#' library('TreeTools')
#' 
#' PathDist(BalancedTree(7), PectinateTree(7))
#' 
#' PathDist(BalancedTree(7), as.phylo(0:2, 7))
#' PathDist(as.phylo(0:2, 7), PectinateTree(7))
#'
#' PathDist(list(bal = BalancedTree(7), pec = PectinateTree(7)),
#'         as.phylo(0:2, 7))
#'
#' CompareAll(as.phylo(30:33, 8), PathDist)
#'   
#' @template MRS
#' @family tree distances
#' @importFrom phangorn path.dist
#' @export
PathDist <- function (tree1, tree2 = tree1) {
  if (inherits(tree1, 'phylo')) {
    tree1 <- Postorder(tree1)
  } else {
    tree1 <- structure(lapply(tree1, Postorder), class = 'multiPhylo')
  }
  
  if (inherits(tree2, 'phylo')) {
    tree2 <- Postorder(tree2)
  } else {
    tree2 <- structure(lapply(tree2, Postorder), class = 'multiPhylo')
  }
  path.dist(tree1, tree2)
}
