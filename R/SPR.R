#' @describeIn SPRWarning for SPR rearrangements
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
  nTips <- tree$Nnode + 1
  if (nTips < 3) return (tree)
  edge   <- tree$edge
  parent <- edge[, 1]
  child  <- edge[, 2]
  nEdge <- length(parent)
  if (nTips == 3) return (ape::root(tree, SampleOne(child[parent==max(parent)], len=2L)))
  
  # Pick an edge at random
  allEdges <- seq_len(nEdge - 1L) + 1L # Only include one root edge
  not1 <- !logical(nEdge)
  not1[1] <- FALSE
  if (is.null(edgeToBreak)) {
    edgeToBreak <- SampleOne(allEdges, len=nEdge - 1L)
  } else {
    if (edgeToBreak > nEdge) return(tree, SPRWarning("edgeToBreak > nEdge"))
    if (edgeToBreak < 1) return(tree, SPRWarning("edgeToBreak < 1"))
    if (edgeToBreak == 1) edgeToBreak <- which(parent == parent[1])[-1] # Use other side of root
  }
  brokenEdge <- seq_along(parent) == edgeToBreak
  brokenEdge.parentNode <- parent[edgeToBreak]
  brokenEdge.childNode  <-  child[edgeToBreak]
  
  if (!is.null(mergeEdge)) { # Quick sanity checks
    if (mergeEdge > nEdge) return(SPRWarning(tree, "mergeEdge value > number of edges"))
    if (length(mergeEdge) !=  1) 
        return(SPRWarning(tree, paste0("mergeEdge value ", paste(mergeEdge, collapse='|'),  
               " invalid; must be NULL or a vector of length 1\n")))
  }  
  
  edgesCutAdrift <- DescendantEdges(edgeToBreak, parent, child)
  edgesRemaining <- !edgesCutAdrift & !brokenEdge
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
  
  if (is.null(mergeEdge)) {
    mergeEdge <- which(!nearBrokenEdge & !edgesOnAdriftSegment & not1)
    nCandidates <- length(mergeEdge)
    if (nCandidates > 1) mergeEdge <- SampleOne(mergeEdge, len=nCandidates)
  }
  
  adriftReconnectionEdge <- warning("mergeEdge[whichAdrift]")
  rootedReconnectionEdge <- warning("mergeEdge[!whichAdrift]")
1

  if(nearBrokenEdge[mergeEdge]) return(SPRWarning(tree, "Selected mergeEdge will not change tree topology."))
  
  if (!nearBrokenEdge[adriftReconnectionEdge]) {
    edgesToInvert <- EdgeAncestry(adriftReconnectionEdge, parent, child, stopAt = edgeToBreak) & !brokenEdge
    #### which(edgesToInvert)
    if (any(edgesToInvert)) {
      tmp <- parent[edgesToInvert]
      parent[edgesToInvert] <- child[edgesToInvert]
      child[edgesToInvert] <- tmp
    }
    reconnectionSideEdges <- edgesToInvert
    reconnectionSideEdges[adriftReconnectionEdge] <- TRUE
    
    repurposedDaughterEdge <- brokenEdgeDaughters & reconnectionSideEdges
    spareDaughterEdge      <- brokenEdgeDaughters & !reconnectionSideEdges
    ######Assert(identical(sum(repurposedDaughterEdge), sum(spareDaughterEdge), 1))
    #### which(repurposedDaughterEdge)
    #### which(spareDaughterEdge)
    child[repurposedDaughterEdge] <- child[spareDaughterEdge]
    child[spareDaughterEdge] <- parent[adriftReconnectionEdge]
    ######Assert(parent[spareDaughterEdge] == brokenEdge.childNode)
    parent[adriftReconnectionEdge] <- brokenEdge.childNode
  }
  if (!nearBrokenEdge[rootedReconnectionEdge]) {
    if (breakingRootEdge) {
      parent[brokenRootDaughters] <- brokenEdge.parentNode
      spareNode <- child[brokenEdgeSister]
      child [brokenEdgeSister] <- child[rootedReconnectionEdge]
      parent[c(edgeToBreak, brokenEdgeSister)] <- spareNode
      child[rootedReconnectionEdge] <- spareNode
    } else {
      parent[brokenEdgeSister] <- parent[brokenEdgeParent]
      parent[brokenEdgeParent] <- parent[rootedReconnectionEdge]
      parent[rootedReconnectionEdge] <- brokenEdge.parentNode
    }
  }
  
  ######Assert(identical(unique(table(parent)), 2L))
  ######Assert(identical(unique(table(child)),  1L))
  ####   matrix(c(parent, child), ncol=2)
  
  tree$edge <- OrderEdgesNumberNodes(parent, child, nTips, nEdge)
  tree
}

#' Rooted SPR 
#' @describeIn SPR Perform \acronym{SPR} rearrangement, retaining position of root
#' @importFrom ape root
#' @export
RootedSPR <- function(tree, edgeToBreak = NULL, mergeEdge = NULL) {
  nTips <- tree$Nnode + 1
  if (nTips < 4) return (SPRWarning(tree, 'Fewer than 4 tips'))
  edge   <- tree$edge
  parent <- edge[, 1]
  child  <- edge[, 2]
  rootNode <- parent[1]
  rootEdges <- parent == rootNode
  nEdge <- length(parent)
  rightTree <- DescendantEdges(1, parent, child)
  selectableEdges <- !rootEdges
  if (sum( rightTree) < 4) selectableEdges[ rightTree] <- FALSE
  if (sum(!rightTree) < 4) selectableEdges[!rightTree] <- FALSE
  if (!any(selectableEdges)) return(SPRWarning(tree, 'No opportunity to rearrange tree due to root position'))

  if (is.null(edgeToBreak)) {
    edgeToBreak <- SampleOne(which(selectableEdges), len=nEdge - 2L) # Pick an edge at random
  } else {
    if (edgeToBreak > nEdge) return(SPRWarning(tree, "edgeToBreak > nEdge"))
    if (edgeToBreak < 1) return(SPRWarning(tree, "edgeToBreak < 1"))
    if (rootEdges[edgeToBreak]) return(SPRWarning(tree, "RootedSPR cannot break root edge; try SPR"))
  }
  repeat {
    edgeInRight <- rightTree[edgeToBreak]
    subtreeWithRoot <- if (edgeInRight) rightTree else !rightTree
    subtreeEdges <- !rootEdges & subtreeWithRoot
    if (sum(edgesCutAdrift <- DescendantEdges(edgeToBreak, parent, child)) > 2) break;
    if (sum(subtreeEdges, -edgesCutAdrift) > 2) break; # the edge itself, and somewheres else
    # TODO check that all expected selections are valid
    selectableEdges[edgeToBreak] <- FALSE
    Assert(any(selectableEdges))
    edgeToBreak <- SampleOne(selectableEdges, len=nEdge - 2L)
  }
  brokenEdge <- seq_along(parent) == edgeToBreak
  brokenEdge.parentNode <- parent[edgeToBreak]
  brokenEdge.childNode  <-  child[edgeToBreak]
  
  edgesRemaining <- !edgesCutAdrift & subtreeEdges
  edgesOnAdriftSegment <- edgesCutAdrift | brokenEdge
  
  if (!is.null(mergeEdge)) { # Quick sanity checks
    if (any(mergeEdge > nEdge)) return(SPRWarning(tree, "mergeEdge value > number of edges"))
    if (length(mergeEdge) > 2 || length(mergeEdge) == 0) 
        return(SPRWarning(tree, paste0("mergeEdge value ", paste(mergeEdge, collapse='|'),  
               " invalid; must be NULL or a vector of length 1 or 2\n  ")))
    if (length(mergeEdge) == 2 && mergeEdge[1] == mergeEdge[2]) 
      return(SPRWarning(tree, "mergeEdge values must differ"))
    if (!all(subtreeWithRoot[mergeEdge])) return(SPRWarning(tree, paste("mergeEdge", 
          mergeEdge[1], mergeEdge[2], "not on same side of root as edgeToBreak", edgeToBreak)))
  }  
  
  brokenEdgeParent <- child  == brokenEdge.parentNode
  brokenEdgeSister <- parent == brokenEdge.parentNode & !brokenEdge
  
  brokenEdgeDaughters <- parent == brokenEdge.childNode
  nearBrokenEdge <- brokenEdgeSister | brokenEdgeParent | brokenEdgeDaughters | brokenEdge
  Assert(any(brokenEdgeParent))
  
  if (is.null(mergeEdge)) {
    mergeEdge  <- which(subtreeEdges & !nearBrokenEdge)
    nCandidates <- length(mergeEdge)
    if (nCandidates > 1) mergeEdge <- SampleOne(mergeEdge, len=nCandidates)
  }
  if (length(mergeEdge) == 1) {
    if (edgesOnAdriftSegment[mergeEdge]) {
      adriftReconnectionEdge <- mergeEdge
      if (nearBrokenEdge[mergeEdge]) {
        samplable <- which(subtreeEdges & !edgesOnAdriftSegment & !nearBrokenEdge)
      } else {
        samplable <- which(subtreeEdges & !edgesOnAdriftSegment)
        Assert(length(samplable) > 0)
      }
      nSamplable <- length(samplable)
      if (nSamplable == 0) return(SPRWarning(tree, "No reconnection site would modify the tree; check mergeEdge"))
      rootedReconnectionEdge <- if (nSamplable == 1) samplable else SampleOne(samplable, len=nSamplable)
      #### cat(" - Selected rooted Reconnection Edge: ", rootedReconnectionEdge, "\n")  #### DEBUGGING AID
    } else {
      rootedReconnectionEdge <- mergeEdge
      if (nearBrokenEdge[mergeEdge]) {
        samplable <- which(subtreeEdges & edgesOnAdriftSegment & !nearBrokenEdge)
      } else {
        samplable <- which(subtreeEdges & edgesOnAdriftSegment)
      }
      nSamplable <- length(samplable)
      if (nSamplable == 0) return(SPRWarning(tree, "No reconnection site would modify the tree; check mergeEdge"))
      adriftReconnectionEdge <- if (nSamplable == 1) samplable else SampleOne(samplable, len=nSamplable)
      #### cat(" - Selected adrift Reconnection Edge: ", adriftReconnectionEdge, "\n") #### DEBUGGING AID
    }
  } else {
    whichAdrift <- edgesOnAdriftSegment[mergeEdge]
    if (sum(whichAdrift) != 1) return(SPRWarning(tree, paste("Invalid edges selected to merge:",
            mergeEdge[1], mergeEdge[2])))
    adriftReconnectionEdge <- mergeEdge[whichAdrift]
    rootedReconnectionEdge <- mergeEdge[!whichAdrift]
  }
  if(nearBrokenEdge[rootedReconnectionEdge] && nearBrokenEdge[adriftReconnectionEdge]) 
    return(SPRWarning(tree, "Selected mergeEdge will not change tree topology."))
  #### edgelabels(edge = edgeToBreak, bg='orange', cex=1.8)  #### DEBUGGING AID
  #### edgelabels(edge=adriftReconnectionEdge, bg='cyan')    #### DEBUGGING AID
  #### edgelabels(edge=rootedReconnectionEdge, bg='magenta') #### DEBUGGING AID
  
  Assert(edgesOnAdriftSegment[adriftReconnectionEdge])
  Assert(!edgesOnAdriftSegment[rootedReconnectionEdge])
  
  if (!nearBrokenEdge[adriftReconnectionEdge]) {
    edgesToInvert <- EdgeAncestry(adriftReconnectionEdge, parent, child, stopAt = edgeToBreak) & !brokenEdge
    if (any(edgesToInvert)) {
      tmp <- parent[edgesToInvert]
      parent[edgesToInvert] <- child[edgesToInvert]
      child[edgesToInvert] <- tmp
    }
    reconnectionSideEdges <- edgesToInvert
    reconnectionSideEdges[adriftReconnectionEdge] <- TRUE
    
    repurposedDaughterEdge <- brokenEdgeDaughters & reconnectionSideEdges
    spareDaughterEdge      <- brokenEdgeDaughters & !reconnectionSideEdges
    Assert(identical(sum(repurposedDaughterEdge), sum(spareDaughterEdge), 1))
    #### which(repurposedDaughterEdge)
    #### which(spareDaughterEdge)
    child[repurposedDaughterEdge] <- child[spareDaughterEdge]
    child[spareDaughterEdge] <- parent[adriftReconnectionEdge]
    Assert(parent[spareDaughterEdge] == brokenEdge.childNode)
    parent[adriftReconnectionEdge] <- child[edgeToBreak]
  }
  if (!nearBrokenEdge[rootedReconnectionEdge]) {
    parent[brokenEdgeSister] <- parent[brokenEdgeParent]
    parent[brokenEdgeParent] <- parent[rootedReconnectionEdge]
    parent[rootedReconnectionEdge] <- brokenEdge.parentNode
  }
  
  Assert(identical(unique(table(parent)), 2L))
  Assert(identical(unique(table(child)),  1L))
  ####   matrix(c(parent, child), ncol=2)
  
  retTree <- tree
  retTree$edge <- OrderEdgesNumberNodes(parent, child, nTips, nEdge)
  retTree
}
