#' NNI
#'
#' Nearest Neighbour Interchange
#' 
#' Performs a single iteration of the nearest-neighbour interchange algorithm.
#' Based on the corresponding \code{phangorn} function, but re-coded to improve speed.
#' 
#' Branch lengths are not supported.
#' 
#' @template treeParam
#' @template edgeToBreakParam
#' 
#' @return Returns a tree with class \code{phylo}.
#'
#' @references
#' The algorithm is summarized in
#' Felsenstein, J. 2004. \cite{Inferring Phylogenies.} Sinauer Associates, Sunderland, Massachusetts.
#' 
#' @author Martin R. Smith
#' 
#' @examples
#' tree <- ape:::rtree(20, br=NULL)
#' NNI(tree)
#'
#' @export
NNI <- function (tree, edgeToBreak=NULL) {
  edge    <- tree$edge
  parent  <- edge[, 1]
  child   <- edge[, 2]
  nTips  <- (length(parent) / 2L) + 1L
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

#' Rooted NNI 
#' @describeIn NNI Perform \acronym{NNI} rearrangement, retaining position of root
#' @export
RootedNNI <- function (tree, edgeToBreak = NULL) {
  edge    <- tree$edge
  parent  <- edge[, 1]
  child   <- edge[, 2]
  nTips  <- (length(parent) / 2L) + 1L
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