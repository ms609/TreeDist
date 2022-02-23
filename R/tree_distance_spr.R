#' Approximate Subtree Prune and Regraft distance
#' 
#' Approximate the Subtree Prune and Regraft (SPR) distance.
#' 
#' `SPRDist()` is a wrapper for the function 
#' \code{\link[phangorn:treedist]{SPR.dist()}} in the phangorn package.
#' It pre-processes trees to ensure that their internal representation does
#' not cause the `SPR.dist()` function to crash R, and allows an improved
#' (but slower) symmetric heuristic.
#' 
#' Note that the phangorn implementation calculates a lower bound on the SPR,
#' using the method of \insertCite{deOliveira2008;textual}{TreeDist}.
#' Other approximations are available
#' \insertCite{@e.g. @Goloboff2008SPR, @Whidden2018}{TreeDist}.
#' 
#' @template tree12ListParams
#' @param symmetric Ignored (redundant after fix of
#' [phangorn#97](https://github.com/KlausVigo/phangorn/issues/97)).
#' 
#' @return `SPRDist()` returns a vector or distance matrix of distances 
#' between trees.
#' 
#' @references \insertAllCited{}
#' 
#' @examples
#' library('TreeTools', quietly = TRUE, warn.conflicts = FALSE)
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
#' @template MRS
#'   
#' @seealso More sophisticated calculation with [\pkg{TBRDist}](
#' https://ms609.github.io/TBRDist/reference/TreeRearrangementDistances.html)
#' functions `USPRDist()` and `ReplugDist()`.
#' @family tree distances
#' @importFrom phangorn SPR.dist
#' @importFrom TreeTools Postorder
#' @export
SPRDist <- function(tree1, tree2 = NULL, symmetric) {
  if (inherits(tree1, 'phylo')) {
    tree1 <- Postorder(tree1)
  } else {
    if (inherits(tree2, 'multiPhylo')) {
      return(vapply(tree2, SPRDist, double(length(tree1)), tree1))
    }
    tree1 <- structure(lapply(tree1, Postorder), class = 'multiPhylo')
  }
  
  if (inherits(tree2, 'phylo')) {
    tree2 <- Postorder(tree2)
  } else if (!is.null(tree2)) {
    tree2 <- structure(lapply(tree2, Postorder), class = 'multiPhylo')
  }
  
  SPR.dist(tree1, tree2)
}
