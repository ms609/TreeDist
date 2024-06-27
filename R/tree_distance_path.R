#' Path distance
#' 
#' Calculate the path distance between rooted or unrooted trees.
#' 
#' This function is a faster alternative to the function 
#' \code{\link[phangorn:treedist]{path.dist()}} in the phangorn package,
#' which can crash if the internal representation of trees does not conform to
#' certain (unspecified) expectations, and which treats all trees as unrooted.
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
#' PathDist(as.phylo(30:33, 8))
#'  
#' @references \insertAllCited{}
#' 
#' @template MRS
#' @family tree distances
#' @importFrom TreeTools Postorder
#' @export
PathDist <- function(tree1, tree2 = NULL) {
  if (inherits(tree1, "phylo")) {
    if (inherits(tree2, "phylo")) {
      .PathDist11(tree1, tree2)
    } else {
      .PathDist1Many(tree1, tree2)
    }
  } else if (is.null(tree2)) {
    .PathDistManySelf(tree1)
  } else if (inherits(tree2, "phylo")) {
    .PathDist1Many(tree2, tree1)
  } else {
    .PathDistManyMany(tree1, tree2)
  }
}

.EuclideanDistance <- function(x) sqrt(sum(x * x))

.PathDist11 <- function(tree1, tree2) {
  .EuclideanDistance(PathVector(tree1) - PathVector(RenumberTips(tree2, tree1)))
}

.PathDist1Many <- function(tree1, treeMany) {
  v1 <- PathVector(tree1)
  apply(v1 - vapply(RenumberTips(treeMany, tree1), PathVector, v1), 2,
        .EuclideanDistance)
}

.PathDistManyMany <- function(trees1, trees2) {
  nTip <- NTip(trees1[[1]])
  v1 <- vapply(RenumberTips(trees1, trees1), PathVector,
               integer(nTip * (nTip - 1) / 2))
  v2 <- vapply(RenumberTips(trees2, trees1), PathVector,
               integer(nTip * (nTip - 1) / 2))
  vec_diff_euclidean(v1, v2)
}

.PathDistManySelf <- function(trees) {
  nTip <- NTip(trees[[1]])
  v1 <- vapply(RenumberTips(trees, trees), PathVector,
               integer(nTip * (nTip - 1) / 2))
  
  nTree <- length(trees)
  
  ret <- structure(pair_diff_euclidean(v1),
                   Size = nTree, Diag = FALSE, Upper = FALSE,
                   class = "dist")
  
  # Return:
  ret
}
