#' Approximate Subtree Prune and Regraft distance
#' 
#' Approximate the Subtree Prune and Regraft (SPR) distance.
#' 
#' This function is a wrapper for the function 
#' \code{\link[phangorn:treedist]{SPR.dist()}} in the phangorn package.
#' It pre-processes trees to ensure that their internal representation does
#' not cause the `SPR.dist()` function to crash R.
#' 
#' @template tree12ListParams
#' @param symmetric Logical specifying whether to produce a better heuristic
#' by calculating the minimum of `SPRDist(t1, t2)` and `SPRDist(t2, t1)`,
#' which are not guaranteed to be equal due to the heuristic nature of the 
#' approximation (see 
#' [phangorn#97](https://github.com/KlausVigo/phangorn/issues/97)). Set to
#' `FALSE` for the faster approximation, as implemented in 'phangorn'.
#' 
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
SPRDist <- function (tree1, tree2 = NULL, symmetric = TRUE) {
  if (inherits(tree1, 'phylo')) {
    tree1 <- Postorder(tree1)
  } else {
    if (inherits(tree2, 'multiPhylo')) {
      return(vapply(tree2, SPRDist, double(length(tree1)), tree1, symmetric))
    }
    tree1 <- structure(lapply(tree1, Postorder), class = 'multiPhylo')
  }
  
  if (inherits(tree2, 'phylo')) {
    tree2 <- Postorder(tree2)
  } else if (!is.null(tree2)) {
    tree2 <- structure(lapply(tree2, Postorder), class = 'multiPhylo')
  }
  
  if (symmetric) {
    if (is.null(tree2)) {
      pmin(SPR.dist(tree1), SPR.dist(rev(tree1)))
    } else {
      if (inherits(tree1, 'phylo') && inherits(tree2, 'multiPhylo')) {
        backwards <- vapply(tree2, SPR.dist, 0, tree1)
      } else if (inherits(tree2, 'phylo') && inherits(tree1, 'multiPhylo')) {
        backwards <- vapply(tree1, SPR.dist, 0, tree2)
      } else {
        backwards <- SPR.dist(tree2, tree1)
      }
      pmin(SPR.dist(tree1, tree2), backwards)
    }
  } else {
    SPR.dist(tree1, tree2)
  }
}
