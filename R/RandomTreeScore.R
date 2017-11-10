#' Parsimony score of random postorder tree
#' 
#' @param nTip number of tips (minimum 3)
#' @template morphyObjParam
#'
#' @return the parsimony score of a random tree, for the given Morphy object.
#'
#' @export
RandomTreeScore <- function (nTip, morphyObj) {  
  if (nTip < 3) {
    warning("nTip < 3 not implemented, as there's only one unrooted topology.")
    return (0)
  }
  # Return:
  .Call('RANDOM_TREE_SCORE', as.integer(nTip), morphyObj)
}

#' Random postorder tree
#' 
#' @param nTip number of tips (minimum 3)
#'
#' @return A list with three elements, each a vector of integers, respectively containing:
#'         - The parent of each tip and node, in order
#'         - The left child of each node
#'         - The right child of each node.
#'
#' @export
RandomMorphyTree <- function (nTip) {  
  if (nTip < 2) {
    stop("nTip < 2 not implemented: a tip is not a tree.")
  }
  # Return:
  .Call('RANDOM_TREE', as.integer(nTip))
}

#' @importFrom graphics plot
plot.morphyTree <- function (morphyTree) {
  parentOf <- morphyTree[[1]]
  left <- morphyTree[[2]]
  right <- morphyTree[[3]]
  nTip <- length(left) + 1L
  
  edge <- matrix(c(rep(seq(nTip, len=nTip - 1L), 2), right, left), ncol=2) + 1L
  tree <- list(edge=edge, Nnode=nTip - 1L, tip.label=seq_len(nTip) - 1L)
  class(tree) <- 'phylo'
  plot(tree)
}
