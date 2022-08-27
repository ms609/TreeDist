#' Path distance
#' 
#' Calculate the path distance between rooted or unrooted trees.
#' 
#' This function is a wrapper for the function 
#' \code{\link[phangorn:treedist]{path.dist()}} in the phangorn package.
#' It pre-processes trees to ensure that their internal representation does
#' not cause the `path.dist()` function to crash R.
#' 
#' The path distance is calculated by tabulating the cladistic difference (=
#' topological distance) between each pair of tips in each tree.
#' A precursor to the path distance \insertCite{Farris1969}{TreeDist}
#' took the mean squared 
#' difference between the elements of each tree's tabulation (Farris, 1973);
#' the method used here is that proposed by
#' \insertCite{Steel1993;textual}{TreeDist}, which takes the square root of this
#' sum.
#' Other precursor measures are described in 
#' \insertCite{Williams1971;textual}{TreeDist} and
#' \insertCite{Phipps1971;textual}{TreeDist}.
#' 
#' If a root node is present, trees are treated as rooted.
#' To avoid counting the root edge twice, use `UnrootTree(tree)` before
#' calculating the path distance.
#' 
#' Use of the path distance is discouraged as it emphasizes 
#' shallow relationships at the expense of deeper (and arguably more
#' fundamental) relationships \insertCite{Farris1973}{TreeDist}.
#' 
#' @template tree12ListParams
#' 
#' @return `PathDist()` returns a vector or distance matrix of distances
#' between trees.
#' 
#' @examples
#' library("TreeTools")
#' 
#' # Treating the two edges to the root node as distinct
#' PathDist(BalancedTree(7), PectinateTree(7))
#' 
#' # Counting those two edges once
#' PathDist(UnrootTree(BalancedTree(7)), UnrootTree(PectinateTree(7)))
#' 
#' PathDist(BalancedTree(7), as.phylo(0:2, 7))
#' PathDist(as.phylo(0:2, 7), PectinateTree(7))
#'
#' PathDist(list(bal = BalancedTree(7), pec = PectinateTree(7)),
#'         as.phylo(0:2, 7))
#'
#' CompareAll(as.phylo(30:33, 8), PathDist)
#'  
#' @references \insertAllCited{}
#' 
#' @template MRS
#' @family tree distances
#' @importFrom phangorn path.dist
#' @importFrom TreeTools Postorder
#' @export
PathDist <- function(tree1, tree2 = NULL) {
  if (inherits(tree1, "phylo")) {
    tree1 <- Postorder(tree1)
  } else {
    tree1 <- structure(lapply(tree1, Postorder), class = "multiPhylo")
  }
  
  if (inherits(tree2, "phylo")) {
    tree2 <- Postorder(tree2)
  } else if (!is.null(tree2)) {
    tree2 <- structure(lapply(tree2, Postorder), class = "multiPhylo")
  }
  path.dist(tree1, tree2)
}
