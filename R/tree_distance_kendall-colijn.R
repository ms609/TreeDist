#' Kendall-Colijn distance between unrooted topologies
#' 
#' The Kendall-Colijn distance is related to the path difference. 
#' 
#' It works by measuring, for each pair of tips, the distance from the most recent
#' common ancestor of those tips and the root node.  For a given tree, this 
#' produces a vector of values recording the distance-from-the-root of each
#' most recent common ancestor of each pair of tips.
#' 
#' Two trees are compared by taking the Euclidian distance between the
#' respective vectors.  This is calculated by taking the square root of the sum 
#' of the squares of the differences between the vectors.
#' 
#' This metric emphasizes the position of the root; the path difference 
#' instead measures the distance of the last common ancestor of each pair
#' of tips from the tips themselves, i.e. the length of the path from one 
#' tip to another.
#' 
#' @seealso [treespace::treeDist](https://cran.r-project.org/web/packages/treespace/vignettes/introduction.html),
#' a more sophisticated, if more cumbersome, implementation that supports 
#' lambda > 0, i.e. use of edge lengths in tree comparison.
#' 
#' @author Martin R. Smith
#' @references \insertRef{Kendall2016}{TreeDist}
#' @export
KendallColijn <- function (tree1, tree2) 
{
  if (!setequal(tree1$tip.label, tree2$tip.label)) 
    stop("Tree tips must bear identical labels")
  
  # Return:
  sqrt(sum((KCVector(tree1) - KCVector(tree2)) ^ 2))
}

#' @describeIn KendallColijn Creates vectors that characterise a rooted tree
#' @importFrom TreeSearch AllAncestors Preorder
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

