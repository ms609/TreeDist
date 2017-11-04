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
#' @return Returns a tree with class \code{phylo} (if \code{returnAll = FALSE}) or 
#'         a set of trees, with class \code{multiPhylo} (if \code{returnAll = TRUE}).
#'
#' @references
#' The algorithm is summarized in
#' Felsenstein, J. 2004. \cite{Inferring Phylogenies.} Sinauer Associates, Sunderland, Massachusetts.
#' 
#' @author Martin R. Smith
#' 
#' @examples
#' tree <- ape::rtree(20, br=NULL)
#' NNI(tree)
#' NNI(tree, all=TRUE)
#'
#' @export
NNI <- function (tree, edgeToBreak=NULL) {
  edge    <- tree$edge
  parent  <- edge[, 1]
  child   <- edge[, 2]
  nEdge   <- length(parent)
  nTips   <- (nEdge / 2L) + 1L
  rootNode <- nTips + 1L
  
  if (!any(samplable)) stop("Not enough edges to allow NNI rearrangement")
  
  if (is.null(edgeToBreak)) { 
    edgeToBreak <- SampleOne(which(samplable))
  } else if (edgeToBreak == -1) {
    # newEdges <- vapply(which(samplable), DoubleNNI, parent=parent, child=child, list(matrix(0L, nEdge, 2), matrix(0L, nEdge, 2)))
    newEdges <- unlist(lapply(which(samplable), DoubleNNI, parent=parent, child=child), recursive=FALSE) # Quicker than vapply, surprisingly
    newTrees <- lapply(newEdges, function (edges) {tree$edge <- edges; tree}) # Quicker than vapply, surprisingly
    return(newTrees)
  } else if (!samplable[edgeToBreak]) {
    stop("edgeToBreak must be an internal edge")
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
  return(tree)
}

#' Double NNI
#' 
#' Returns the edge parameter of the two trees consistent with the speficied \acronym{NNI} rearrangement
#'
#' @param parent tree$edge[, 1]
#' @param child  tree$edge[, 2]
#' @param edgeToBreak index of edge to break
#'
#' @return the \code{tree$edge} parameter of the two trees consistent with the specified rearrangement
#'
#' @keywords internal
#' @author Martin R. Smith
#'
DoubleNNI <- function (parent, child, edgeToBreak) {
  end1   <- parent[edgeToBreak]
  end2   <- child[edgeToBreak]
  ind1   <- which(parent == end1)
  ind1   <- ind1[ind1 != edgeToBreak][1]
  ind2.3 <- which(parent == end2)
  ind2   <- ind2.3[1]
  ind3   <- ind2.3[2]

  newInd <- c(ind2, ind1)
  oldInd <- c(ind1, ind2)
  child2 <- child
  childSwap <- child[newInd]
  child2[oldInd] <- childSwap
  
  newInd <- c(ind3, ind1)
  oldInd <- c(ind1, ind3)
  childSwap <- child[newInd]
  child[oldInd] <- childSwap
  
  nEdge <- length(parent)
  return(list(RenumberTree(parent, child, nEdge), RenumberTree(parent, child2, nEdge)))
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
  } else if (edgeToBreak == -1) {
    newEdges <- unlist(lapply(which(samplable), DoubleNNI, parent=parent, child=child), recursive=FALSE) # Quicker than vapply, surprisingly
    newTrees <- lapply(newEdges, function (edges) {tree$edge <- edges; tree}) # Quicker than vapply, surprisingly
    return(newTrees)   
  } else if (!samplable[edgeToBreak]) {
    stop("edgeToBreak cannot include a tip or the root node")
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