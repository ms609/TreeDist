#' Kendall-Colijn distance between unrooted topologies
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
  
  tip2 <- matrix(tipOrder, nrow=nTip, ncol=nTip)
  ancestors <- AllAncestors(parent, child)
  
  mrca <- mapply(function(anc1, anc2) max(intersect(anc1, anc2)), 
                 ancestors[rep(tipOrder, nTip - seq_len(nTip))],
                 ancestors[tip2[lower.tri(tip2)]])
  
  rootDist <- vapply(ancestors, length, integer(1))
  rootDist[mrca]
}
