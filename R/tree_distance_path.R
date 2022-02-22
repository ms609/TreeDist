#' Path distance
#' 
#' Calculate the path distance between rooted or unrooted trees.
#' 
#' This function is an alternative to the function 
#' \code{\link[phangorn:treedist]{path.dist()}} in the phangorn package,
#' which can crash if the internal representation of trees does not conform to
#' certain (unspecified) expectations.
#' 
#' The path distance is calculated by tabulating the cladistic difference (=
#' topological distance) between each pair of tips in each tree.
#' A precursor to the path distance (Farris, 1969) took the mean squared 
#' difference between the elements of each tree's tabulation (Farris, 1973);
#' the method used here is that proposed by Steel & Penny (1993), which takes
#' the square root of this sum.  Other precursor measures are described in 
#' Williams and Clifford (1971) and Phipps (1971).
#' 
#' If a root node is present, trees are treated as rooted.
#' To avoid counting the root edge twice, use `UnrootTree(tree)` before
#' calculating the path distance.
#' 
#' Use of the path distance is discouraged as it emphasizes 
#' shallow relationships at the expense of deeper (and arguably more
#' fundamental) relationships (Farris, 1973).
#' 
#' @template tree12ListParams
#' 
#' @return `PathDist()` returns a vector or distance matrix of distances
#' between trees.
#' 
#' @examples
#' library('TreeTools')
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
#' @references 
#' \insertRef{Farris1969}{TreeDist}
#' 
#' \insertRef{Farris1973}{TreeDist}
#' 
#' \insertRef{Phipps1971}{TreeDist}
#' 
#' \insertRef{Steel1993}{TreeDist}
#' 
#' \insertRef{Williams1971}{TreeDist}
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
  apply(v2, 2, function(X) apply(X - v1, 2, .EuclideanDistance))
}

.PathDistManySelf <- function(trees) {
  nTip <- NTip(trees[[1]])
  v1 <- vapply(RenumberTips(trees, trees), PathVector,
               integer(nTip * (nTip - 1) / 2))
  
  nTree <- length(trees)
  nPair <- nTree * (nTree - 1) / 2
  ret <- structure(numeric(nPair), Size = nTree, class = "dist",
                   Diag = FALSE, Upper = TRUE)
  ptr <- 0L
  for (i in seq_len(nTree - 1L)) {
    X <- v1[, i]
    entries <- seq_len(nTree - i)
    j <- i + entries
    ret[ptr + entries] <- apply(X - v1[, j, drop = FALSE], 2, .EuclideanDistance)
    ptr <- ptr + length(entries)
  }
  
  # Return:
  ret
}
