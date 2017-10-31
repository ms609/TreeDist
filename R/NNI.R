#' @title Tree rearrangement functions
#' 
#' These functions performs a single random \acronym{TBR}, \acronym{SPR} or \acronym{NNI} iteration.
#'
#' Performs a single iteration of the nearest-neigbour interchange, subtree pruning and regrafting,
#' or tree bisection and reconnection algorithms.
#' NNI and SPR are based on the corresponding phangorn functions, but have been re-coded to 
#' improve their speed.
#' 
#' Branch lengths are not supported.
#' 
#' @usage 
#'  NNI(tree, edgeToBreak = NULL)
#'  SPR(tree, edgeToBreak = NULL, mergeEdges = NULL)
#'  TBR(tree, edgeToBreak = NULL, mergeEdges = NULL)
#'
#' @template treeParam
#' @template edgeBreakingParams
#' 
#' @return Returns a tree with class \code{phylo}.
#'
#' @references
#' The algorithms are summarized in
#' Felsenstein, J. 2004. \cite{Inferring Phylogenies.} Sinauer Associates, Sunderland, Massachusetts.
#' 
#' @author Martin R. Smith
#' 
#' @examples
#' library(ape)
#' tree <- ape:::rtree(20, br=NULL)
#' NNI(tree)
#' SPR(tree)
#' TBR(tree)
#' @export
NNI <- function (tree, edgeToBreak=NULL) {
  edge    <- tree$edge
  parent  <- edge[, 1]
  child   <- edge[, 2]
  nTips  <- length(tree$tip.label)
  rootNode <- nTips + 1L
  
  samplable <- child > nTips
  if (is.null(edgeToBreak)) { 
    edgeToBreak <- SampleOne(which(samplable))
  } else {
    if (!samplable[edgeToBreak]) stop("edgeToBreak must be an internal edge")
  }
  if (is.na(edgeToBreak)) stop("Cannot find a valid rearrangement")
  
  end1   <- parent[edgeToBreak]
  end2   <- child[edgeToBreak]
  ind1   <- which(parent == end1)
  ind1   <- ind1[ind1 != edgeToBreak][1]
  ind2   <- which(parent == end2)[sample.int(2L, 1L, useHash=FALSE)]

  newInd <- c(ind2, ind1)
  oldInd <- c(ind1, ind2)
  childSwap <- child[newInd]
  child[oldInd] <- childSwap
  
  tree$edge <- RenumberTree(parent, child)
  tree
}

#' Rearrange a rooted tree
#'
#' This function performs a rearrangement iteration on a tree, retaining the position of the root.
#'
#' A single \acronym{NNI}, \acronym{SPR} or \acronym{TBR} rearrangement is performed, subject to the constraint that 
#' no taxon may be moved to the opposite side of the root node.
#' Branch lengths are not (yet) supported.
#' 
#' @usage
#' RootedNNI(tree)
#' RootedSPR(tree)
#' RootedTBR(tree)
#'
#' @param tree A bifurcating tree of class \code{\link{phylo}}, with all nodes resolved
#' @template edgeBreakingParams
#' 
#' @return This function returns a tree, in \code{phylo} format.
#'
#' @author Martin R. Smith
#' \code{RootedNNI} is abridged from the \pkg{phangorn} function \code{nnin}
#' 
#' @seealso
#' \itemize{
#' \item \code{\link{NNI}}, unrooted \acronym{NNI} and \acronym{SPR}
#' \item \code{\link{TBR}}, unrooted \acronym{TBR}
#' }
#' 
#' @examples{
#'   require('ape')
#'   tree <- read.tree(text='(((a,b),c),(d,(e,f)));')
#'   tree <- root(tree, c('e', 'f'), resolve.root=TRUE)
#'   plot(tree)
#'   dev.new()
#'   plot(RootedNNI(tree))
#'   plot(RootedSPR(tree))
#'   plot(RootedTBR(tree))
#' }
#' 
#'
#' @export
RootedNNI <- function (tree, edgeToBreak = NULL) {
  edge    <- tree$edge
  parent  <- edge[, 1]
  child   <- edge[, 2]
  nTips  <- length(tree$tip.label)
  rootNode <- nTips + 1L
  
  samplable <- parent != rootNode & child > nTips
  if (is.null(edgeToBreak)) { 
    edgeToBreak <- SampleOne(which(samplable))
  } else {
    if (!samplable[edgeToBreak]) stop("edgeToBreak cannot include a tip or the root node")
  }
  if (is.na(edgeToBreak)) stop("Cannot find a valid rearrangement")
  
  end1   <- parent[edgeToBreak]
  end2   <- child[edgeToBreak]
  ind1   <- which(parent == end1)
  ind1   <- ind1[ind1 != edgeToBreak][1]
  ind2   <- which(parent == end2)[sample.int(2L, 1L, useHash=FALSE)]
  
  newInd <- c(ind2, ind1)
  oldInd <- c(ind1, ind2)
  
  child_swap <- child[newInd]
  child[oldInd] <- child_swap
  tree$edge <- RenumberTree(parent, child)
  tree
}