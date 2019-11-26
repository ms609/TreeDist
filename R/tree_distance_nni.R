# An approximation of the NNI distance based on Li et al. 1996

# Ultimately this should be moved to raw-data.
SortingNumber <- function (n) {
  logCeiling <- ceiling(log2(n))
  n * logCeiling - (2L ^ logCeiling) + 1
}

# Ultimately this should be moved to raw-data.
DegenerateDistance <- function (nTip) {
  # Calculate the shortest length of the longest path in a tree.
  # To make this path as short as possible, divide tips into three 
  # balanced trees, joined by a single node that will form part of every 
  # longest path.  One of these subtrees will be filled with >= n/3 nodes
  nodesInFull <- pmax(0, ceiling(log2(nTip / 3)))
  # We want to put a power of two tips in this subtree, such that every node is 
  # equally close to its root
  tipsInFull <- 2 ^ nodesInFull
  # Now the remaining tips must be spread sub-evenly between the remaining 
  # edges from this node.  Picture halving the tips; removing tips from one side
  # until it is a power of two will reduce the number of nodes by one, whilst 
  # at worst (if this brings the other side over a power of two) increasing the 
  # other side by one.
  tipsLeft <- nTip - tipsInFull
  minBackboneNodes <- pmax(0, nodesInFull + ceiling(log2(tipsLeft / 2)) + 1L)
  # The worst-case scenario requires a move for every node not on the backbone:
  nNode <- pmax(0, nTip - 2L)
  
  # Return:
  nNode - minBackboneNodes
}

MaxNNI <- function (nEdge) {
  if (nEdge < 12L) {
    # Exact value of diameter for small trees, calculated by Li et al. 1996.
    c(0L, 0L, 0L, 1L, 3L, 5L, 7L, 10L, 12L, 15L, 18L)[nEdge]
  } else {
    
    nTip <- nEdge + 2L # not 3L??
    # Unclear from Li et al. 1996 whether the DegenerateDistance must be applied
    # twice.  Conservatively assumed so: once to get to degenerate tree, and once
    # to get back to balanced tree.
    SortingNumber(nEdge + 3L) + (2L * DegenerateDistance(nTip))
  }
}

MinNNI <- function (n) n

#' NNI distance
#' 
#' @examples
#'   library('TreeTools')
#'   tree1 <- BalancedTree(7)
#'   tree2 <- PectinateTree(7)
#'   
#'   NNIDist(tree1, tree2)
#' 
#' @template MRS
#' @importFrom ape Nnode.phylo
#' @importFrom TreeTools PostorderEdges NTip RenumberTips
#' @export
NNIDist <- function (tree1, tree2) {
  edge <- tree1$edge
  edge1 <- do.call(cbind, 
                   PostorderEdges(edge[, 1], edge[, 2], dim(edge)[1], 
                                  Nnode.phylo(tree1)))
  tree2 <- RenumberTips(tree2, tree1$tip.label)
  edge <- tree2$edge
  edge2 <- do.call(cbind, 
                   PostorderEdges(edge[, 1], edge[, 2], dim(edge)[1], 
                                  Nnode.phylo(tree2)))
  cpp_nni_distance(edge1, edge2, NTip(tree1))
}