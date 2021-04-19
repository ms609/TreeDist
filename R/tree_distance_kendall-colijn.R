#' Kendall&ndash;Colijn distance
#' 
#' Calculate the Kendall&ndash;Colijn tree distance, a measure related to the 
#' path difference. 
#' 
#' The Kendall&ndash;Colijn distance works by measuring, for each pair of leaves,
#' the distance from the most recent
#' common ancestor of those leaves and the root node.  For a given tree, this 
#' produces a vector of values recording the distance-from-the-root of each
#' most recent common ancestor of each pair of leaves.
#' 
#' Two trees are compared by taking the Euclidian distance between the
#' respective vectors.  This is calculated by taking the square root of the sum 
#' of the squares of the differences between the vectors.
#' 
#' This metric emphasizes the position of the root; the path difference 
#' instead measures the distance of the last common ancestor of each pair
#' of leaves from the leaves themselves, i.e. the length of the path from one 
#' leaf to another.
#' 
#' @template tree12ListParams
#' 
#' @templateVar returns `KendallColijn()` returns
#' @template distReturn
#' 
#' @examples 
#' KendallColijn(TreeTools::BalancedTree(8), TreeTools::PectinateTree(8))
#'
#' set.seed(0)
#' KendallColijn(TreeTools::BalancedTree(8), lapply(rep(8, 3), ape::rtree))
#' KendallColijn(lapply(rep(8, 4), ape::rtree))
#' @template MRS
#' 
#' @seealso [`treespace::treeDist`](https://CRAN.R-project.org/package=treespace/vignettes/introduction.html)
#' is a more sophisticated, if more cumbersome, implementation that supports 
#' lambda > 0, i.e. use of edge lengths in tree comparison.
#' 
#' @references \insertRef{Kendall2016}{TreeDist}
#' @family tree distances
#' @encoding UTF-8
#' @export
KendallColijn <- function (tree1, tree2 = tree1) {
  FunValue <- function (nTip) double(nTip * (nTip - 1L) / 2L)
  
  EuclidianDistance <- function (x) sqrt(sum(x * x))
  
  if (inherits(tree1, 'phylo')) {
    if (inherits(tree2, 'phylo')) {
      if (length(tree1$tip.label) != length(tree2$tip.label) || 
          length(setdiff(tree1$tip.label, tree2$tip.label)) > 0) {
        stop("Leaves must bear identical labels.")
      }
      EuclidianDistance(KCVector(tree1) - KCVector(tree2))
    } else {
      apply(KCVector(tree1) - vapply(tree2, KCVector,
                                     FunValue(length(tree1$tip.label))),
            2L, EuclidianDistance)
    }
  } else {
    if (inherits(tree2, 'phylo')) {
      apply(KCVector(tree2) - vapply(tree1, KCVector,
                                     FunValue(length(tree2$tip.label))),
            2L, EuclidianDistance)
    } else {
      nTip <- length(tree1[[1]]$tip.label)
      vector1 <- vapply(tree1, KCVector, FunValue(length(tree1[[1]]$tip.label)))
      vector2 <- vapply(tree2, KCVector, FunValue(length(tree2[[1]]$tip.label)))
      apply(vector2, 2, function(i) 
        apply(vector1, 2, function (j) 
          EuclidianDistance(i - j)))
    }
  }
}

#' @describeIn KendallColijn Creates a vector that characterises a rooted tree,
#' as described in Kendall & Colijn (2016).
#' @param tree A tree of class \code{\link[ape:read.tree]{phylo}}.
#' @importFrom TreeTools AllAncestors Preorder
#' @export
KCVector <- function (tree) {
  tree <- Preorder(tree)
  edge <- tree$edge
  parent <- edge[, 1L]
  child <- edge[, 2L]
  root <- parent[1]
  nTip <- root - 1L
  tipEdges <- match(seq_len(nTip), child)
  tipAncs <- seq_len(nTip)
  tipOrder <- order(tree$tip.label)
  
  tipI <- rep(tipOrder, nTip - seq_len(nTip))
  tipJ <- matrix(tipOrder, nrow=nTip, ncol=nTip)
  tipJ <- tipJ[lower.tri(tipJ)]
  
  ancestors <- AllAncestors(parent, child)
  
  mrca <- mapply(function(anc1, anc2) max(intersect(anc1, anc2)), 
                 ancestors[tipI], ancestors[tipJ])
  
  rootDist <- vapply(ancestors, length, integer(1))
  rootDist[mrca]
}

