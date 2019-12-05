#' Approximate the NNI distance
#' 
#' An approximation of the NNI distance based on Li et al. (1996).
#' 
#' In brief, this approximation algorithm works by identifying edges in one
#' tree that do not match edges in the second.  Each of these edges must
#' undergo at least one NNI operation in order to reconcile the trees.
#' Edges that match in both trees need never undergo an NNI operation, and 
#' divide each tree into smaller regions.  By 'cutting' matched edges into two,
#' a tree can be divided into a number of regions that solely comprise unmatched
#' edges.
#' 
#' These regions can be viewed as separate trees that need to be reconciled.
#' One way to reconcile these trees is to conduct a series of NNI operations
#' that reduce a tree to a pectinate (caterpillar) tree, then to conduct an 
#' analogue of the mergesort algorithm.  This takes at most _n_ log _n_ + O(_n_)
#' NNI operations, and provides a loose upper bound on the NNI score. 
#' The maximum number of moves for a tree with _n_ tips can be calculated 
#' exactly for small trees; this provides a tighter upper bound, but is 
#' unavailable for _n_ > 12.
#' 
#' \tabular{rccccccccccccc}{
#'   Tips:     \tab 1 \tab 2 \tab 3 \tab 4 \tab 5 \tab 6 \tab 7 \tab 8 \tab 9
#'    \tab 10 \tab 11 \tab 12 \tab 13 \cr
#'   Diameter: \tab 0 \tab 0 \tab 0 \tab 1 \tab 3 \tab 5 \tab 7 \tab 10 \tab 12
#'    \tab 15 \tab 18 \tab 21 \tab ? \cr
#'  }
#' 
#' @param tree1,tree2 Individual trees of class `phylo`
#' 
#' @return `NNIDist` returns, for each pair of trees, a list containing three
#' entries:
#' 
#' - `lower` is a lower bound on the NNI distance, and corresponds
#' to the RF distance between the trees. 
#' 
#' - `tight_upper` is an upper bound on the distance, based on calculated
#' maximum diameters for trees with <12 tips.  _NA_ is returned if trees are
#' too different to employ this approach.
#' 
#' - `loose_upper` is a looser upper bound on the distance, using _n_ log _n_ +
#' O(_n_).
#' 
#' 
#' @references 
#'   \insertRef{Li1996}{TreeDist}
#' 
#' @examples
#'   library('TreeTools')
#'   tree1 <- BalancedTree(7)
#'   tree2 <- PectinateTree(7)
#'   
#'   NNIDist(tree1, tree2)
#'   
#' @template MRS
#' @export
NNIDist <- function (tree1, tree2 = tree1) {
  .TreeDistance(.NNIDistSingle, tree1, tree2)
}

#' @importFrom TreeTools PostorderEdges
.EdgesInPostorder <- function (edge, nNode) {
  do.call(cbind, PostorderEdges(edge[, 1], edge[, 2], dim(edge)[1], nNode))
}

#' @importFrom TreeTools RenumberTips
#' @importFrom ape Nnode.phylo
.NNIDistSingle <- function (tree1, tree2, nTip, ...) {
  tree2 <- RenumberTips(tree2, tree1$tip.label)
  
  edge1 <- .EdgesInPostorder(tree1$edge, Nnode.phylo(tree1))
  edge2 <- .EdgesInPostorder(tree2$edge, Nnode.phylo(tree2))
  
  cpp_nni_distance(edge1, edge2, nTip)
}


