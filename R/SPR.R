#' @describeIn TBRWarning for SPR rearrangements
#' @keywords internal
#' @export
SPRWarning <- function (tree, error) {
  warning ("No SPR operation performed.\n  > ", error)
  return(tree)
}

#' Subtree Pruning and Rearrangement 
#'
#' Perform one \acronym{SPR} rearrangement on a tree
#' 
#' Equivalent to phangorn's kSPR, but faster.
#' Note that rearrangements that only change the position of the root WILL be returned by 
#' \code{SPR}.  If the position of the root is irrelevant (as in Fitch parsimony, for example)
#' then this function will occasionally return a functionally equivalent topology.  
#' \code{RootIrrelevantSPR} will search tree space more efficiently in these cases.
#' Branch lengths are not (yet) supported.
#'
#' @template treeParam
#' @param edgeToBreak the index of an edge to bisect, generated randomly if not specified.
#' @return This function returns a tree in \code{phyDat} format that has undergone one \acronym{SPR} iteration.
#' 
#' @references The \acronym{SPR} algorithm is summarized in
#' Felsenstein, J. 2004. \cite{Inferring Phylogenies.} Sinauer Associates, Sunderland, Massachusetts.
#' 
#' @author Martin R. Smith
#' 
#' @seealso RootedSPR useful when the position of the root node should be retained.
#' @seealso TBR
#' @seealso NNI
#' 
#' @examples{
#' tree <- RandomTree(20, br=NULL)
#' SPR(tree)
#' }
#' @importFrom ape root
#' @export
SPR <- function(tree, edgeToBreak = NULL, mergeEdge = NULL) {
  if (is.null(treeOrder <- attr(tree, 'order')) || treeOrder != 'preorder') tree <- Preorder(tree)
  nTips <- tree$Nnode + 1L
  if (nTips < 3L) return (tree)
  edge   <- tree$edge
  parent <- edge[, 1]
  child  <- edge[, 2]
  nEdge <- length(parent)
  if (nTips == 3) return (ape::root(tree, SampleOne(child[parent==max(parent)], len=2L)))
  
  notDuplicateRoot <- !logical(nEdge)
  rightSide <- DescendantEdges(1, parent, child, nEdge)
  nEdgeRight <- sum(rightSide)
  if (nEdgeRight == 1) {
    notDuplicateRoot[2] <- FALSE
  } else if (nEdgeRight == 3) {
    notDuplicateRoot[4] <- FALSE
  } else {
    notDuplicateRoot[1] <- FALSE
  }
  
  if (is.null(edgeToBreak)) {
    # Pick an edge at random
    edgeToBreak <- SampleOne(which(notDuplicateRoot), len=nEdge - 1L)
  } else {
    if (edgeToBreak > nEdge) return(SPRWarning(tree, "edgeToBreak > nEdge"))
    if (edgeToBreak < 1) return(SPRWarning(tree, "edgeToBreak < 1"))
  }
  brokenEdge <- seq_along(parent) == edgeToBreak
  brokenEdge.parentNode <- parent[edgeToBreak]
  brokenEdge.childNode  <-  child[edgeToBreak]
    
  edgesCutAdrift <- DescendantEdges(edgeToBreak, parent, child, nEdge)
  edgesOnAdriftSegment <- edgesCutAdrift | brokenEdge
  
  brokenEdgeParent <- child == brokenEdge.parentNode
  brokenEdgeSister <- parent == brokenEdge.parentNode & !brokenEdge
  brokenEdgeDaughters <- parent == brokenEdge.childNode
  nearBrokenEdge <- brokenEdge | brokenEdgeSister | brokenEdgeParent | brokenEdgeDaughters
  if (breakingRootEdge <- !any(brokenEdgeParent)) { 
    # Edge to break is the Root Node.
    brokenRootDaughters <- parent == child[brokenEdgeSister]
    nearBrokenEdge <- nearBrokenEdge | brokenRootDaughters
  }
    
  if (!is.null(mergeEdge)) { # Quick sanity checks
    if (mergeEdge > nEdge) return(SPRWarning(tree, "mergeEdge value > number of edges"))
    if (length(mergeEdge) !=  1) 
        return(SPRWarning(tree, paste0("mergeEdge value ", paste(mergeEdge, collapse='|'),  
               " invalid; must be NULL or a vector of length 1\n")))
    if(nearBrokenEdge[mergeEdge]) return(SPRWarning(tree, "Selected mergeEdge will not change tree topology."))
    if(DescendantEdges(edgeToBreak, parent, child, nEdge)[mergeEdge]) stop("mergeEdge is within pruned subtree")
  } else {
    mergeEdge <- which(!nearBrokenEdge & !edgesOnAdriftSegment & notDuplicateRoot)
    nCandidates <- length(mergeEdge)
    #####Assert(nCandidates > 0)
    if (nCandidates > 1) mergeEdge <- SampleOne(mergeEdge, len=nCandidates)
  }
  
  if (breakingRootEdge) {
    parent[brokenRootDaughters] <- brokenEdge.parentNode
    spareNode <- child[brokenEdgeSister]
    child [brokenEdgeSister] <- child[mergeEdge]
    parent[brokenEdge | brokenEdgeSister] <- spareNode
    child[mergeEdge] <- spareNode
  } else {
    parent[brokenEdgeSister] <- parent[brokenEdgeParent]
    parent[brokenEdgeParent] <- parent[mergeEdge]
    parent[mergeEdge] <- brokenEdge.parentNode
  }
  
  #####Assert(identical(unique(table(parent)), 2L))
  #####Assert(identical(unique(table(child)),  1L))
  ####   matrix(c(parent, child), ncol=2)
  
  tree$edge <- RenumberTree(parent, child, nEdge)
  tree
}

#' Rooted SPR 
#' @describeIn SPR Perform \acronym{SPR} rearrangement, retaining position of root
#' @importFrom ape root
#' @export
RootedSPR <- function(tree, edgeToBreak = NULL, mergeEdge = NULL) {
  if (is.null(treeOrder <- attr(tree, 'order')) || treeOrder != 'preorder') tree <- Preorder(tree)
  nTips <- tree$Nnode + 1L
  if (nTips < 3L) return (tree)
  edge   <- tree$edge
  parent <- edge[, 1]
  child  <- edge[, 2]
  nEdge <- length(parent)
  rootNode <- parent[1]
  rootEdges <- parent == rootNode
  if (nTips == 3) return (SPRWarning(tree, "Only 3 tips - nothing to rearrange!"))
  
  breakable <- !logical(nEdge) & !rootEdges
  rightSide <- DescendantEdges(1, parent, child, nEdge)
  leftSide  <- !rightSide
  nEdgeRight <- which(rootEdges)[2] - 1
  nEdgeLeft <- nEdge - nEdgeRight
  if (nEdgeRight < 4) {
    if (nEdgeLeft < 4) return(SPRWarning(tree, "No rearrangement possible with this root position."))

    breakable <- breakable & !rightSide
    rightHalfOfLeftSide <- DescendantEdges(nEdgeRight + 2L, parent, child, nEdge)
     leftHalfOfLeftSide <- leftSide & !rightHalfOfLeftSide & !rootEdges
      if (sum(rightHalfOfLeftSide) == 1) breakable[nEdgeRight + 3] <- FALSE
      if (sum( leftHalfOfLeftSide) == 1) breakable[nEdgeRight + 2] <- FALSE
  } else {
    if (nEdgeLeft < 4) {
      breakable <- breakable & rightSide
    } else {
      rightHalfOfLeftSide <- DescendantEdges(nEdgeRight + 2L , parent, child, nEdge)
       leftHalfOfLeftSide <- leftSide & !rightHalfOfLeftSide & !rootEdges
      if (sum(rightHalfOfLeftSide) == 1) breakable[nEdgeRight + 3] <- FALSE
      if (sum( leftHalfOfLeftSide) == 1) breakable[nEdgeRight + 2] <- FALSE
    }
    rightHalfOfRightSide <- DescendantEdges(2L , parent, child, nEdge)
     leftHalfOfRightSide <- rightSide & !rightHalfOfRightSide & !rootEdges
    if (sum(rightHalfOfRightSide) == 1) breakable[3] <- FALSE
    if (sum( leftHalfOfRightSide) == 1) breakable[2] <- FALSE
  }  
  
  if (is.null(edgeToBreak)) {
    # Pick an edge at random
    edgeToBreak <- SampleOne(which(breakable))
  } else {
    if (!breakable[edgeToBreak]) return(SPRWarning(tree, paste("Nowhere to regraft if pruning on edge", edgeToBreak)))
    if (edgeToBreak > nEdge) return(SPRWarning(tree, "edgeToBreak > nEdge"))
    if (edgeToBreak < 1) return(SPRWarning(tree, "edgeToBreak < 1"))
  }
  brokenEdge <- seq_along(parent) == edgeToBreak
  brokenEdge.parentNode <- parent[edgeToBreak]
  brokenEdge.childNode  <-  child[edgeToBreak]
  
  edgesCutAdrift <- DescendantEdges(edgeToBreak, parent, child, nEdge)
  edgesOnAdriftSegment <- edgesCutAdrift | brokenEdge
  
  brokenEdgeParent <- child == brokenEdge.parentNode
  brokenEdgeSister <- parent == brokenEdge.parentNode & !brokenEdge
  brokenEdgeDaughters <- parent == brokenEdge.childNode
  nearBrokenEdge <- brokenEdge | brokenEdgeSister | brokenEdgeParent | brokenEdgeDaughters
    
  if (!is.null(mergeEdge)) { # Quick sanity checks
    if (mergeEdge > nEdge) return(SPRWarning(tree, "mergeEdge value > number of edges"))
    if (length(mergeEdge) !=  1) 
        return(SPRWarning(tree, paste0("mergeEdge value ", paste(mergeEdge, collapse='|'),  
               " invalid; must be NULL or a vector of length 1\n")))
    if(nearBrokenEdge[mergeEdge]) return(SPRWarning(tree, "Selected mergeEdge will not change tree topology."))
    if(DescendantEdges(edgeToBreak, parent, child, nEdge)[mergeEdge]) stop("mergeEdge is within pruned subtree")
  } else {
    edgesOnThisSide <- if (rightSide[edgeToBreak]) rightSide else leftSide
    mergeEdge <- which(edgesOnThisSide & !nearBrokenEdge & !edgesOnAdriftSegment)
    nCandidates <- length(mergeEdge)
    Assert(nCandidates > 0)
    if (nCandidates > 1) mergeEdge <- SampleOne(mergeEdge, len=nCandidates)
  }
  
  parent[brokenEdgeSister] <- parent[brokenEdgeParent]
  parent[brokenEdgeParent] <- parent[mergeEdge]
  parent[mergeEdge] <- brokenEdge.parentNode
  
  #####Assert(identical(unique(table(parent)), 2L))
  #####Assert(identical(unique(table(child)),  1L))
  ####   matrix(c(parent, child), ncol=2)
  
  tree$edge <- RenumberTree(parent, child, nEdge)
  tree
}