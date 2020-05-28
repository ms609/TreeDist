#' Approximate Subtree Prune and Regraft distance
#' 
#' Approximate the Subtree Prune and Regraft (SPR) distance.
#' 
#' This function is a wrapper for the function 
#' \code{\link[phangorn:treedist]{SPR.dist()}} in the phangorn package.
#' It pre-processes trees to ensure that their internal representation does
#' not cause the `SPR.dist()` function to crash R.
#' 
#' A memory leak is present in phangorn v2.5.5.  To avoid a drain on system
#' resources, install the latest version of phangorn with 
#' `devtools::install_github('KlausVigo/phangorn')`.
#' 
#' @template tree12Params
#' 
#' @return `SPRDist()` returns a vector or distance matrix of distances 
#' between trees.
#' 
#' @examples
#' library('TreeTools')
#' 
#' SPRDist(BalancedTree(7), PectinateTree(7))
#' 
#' SPRDist(BalancedTree(7), as.phylo(0:2, 7))
#' SPRDist(as.phylo(0:2, 7), PectinateTree(7))
#'
#' SPRDist(list(bal = BalancedTree(7), pec = PectinateTree(7)),
#'         as.phylo(0:2, 7))
#'
#' CompareAll(as.phylo(30:33, 8), SPRDist)
#'   
#' @template MRS
#' @family tree distances
#' @importFrom phangorn SPR.dist
#' @importFrom TreeTools Postorder
#' @export
SPRDist <- function (tree1, tree2 = NULL) {
  if (inherits(tree1, 'phylo')) {
    tree1 <- Postorder(tree1)
  } else {
    tree1 <- structure(lapply(tree1, Postorder), class = 'multiPhylo')
  }
  
  if (inherits(tree2, 'phylo')) {
    tree2 <- Postorder(tree2)
  } else if (!is.null(tree2)) {
    tree2 <- structure(lapply(tree2, Postorder), class = 'multiPhylo')
  }
  SPR.dist(tree1, tree2)
}
