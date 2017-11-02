#' TBR Warning
#' Print a warning and return given tree
#'
#' @param tree tree to return
#' @param error error message to report
#'
#' @return the tree specified in tree
#' @examples
#' testFunction <- function (tree) {
#'  return(TBRWarning(tree, 'Message text'))
#' }
#' \dontrun{testFunction(0) # will trigger warning}
#' 
#'
#' @author Martin R. Smith
#' @keywords internal
#' @export
TBRWarning <- function (tree, error) {
  warning ("No TBR operation performed.\n  > ", error)
  return(tree)
}

#' TBR
#' 
#' Tree bisection and reconnection
#'
#' \code{TBR} performs a single random \acronym{TBR} iteration.
#'
#' @param tree A bifurcating tree of class \code{\link{phylo}}, with all nodes resolved;
#' @template edgeToBreakParam
#' @template mergeEdgesParam
#' 
#' @details Branch lengths are not (yet) supported.
#' 
#' @return This function returns a tree in \code{phyDat} format that has undergone one \acronym{TBR} iteration.
#' @references The \acronym{TBR} algorithm is summarized in
#' Felsenstein, J. 2004. \cite{Inferring Phylogenies.} Sinauer Associates, Sunderland, Massachusetts.
#' 
#' 
#' @author Martin R. Smith
#' 
#' @seealso RootedTBR useful when the position of the root node should be retained.
#' 
#' @examples{
#' library('ape')
#' tree <- rtree(20, br=NULL)
#' TBR(tree)
#' }
#' @importFrom ape root
#' @export
TBR <- function(tree, edgeToBreak = NULL, mergeEdges = NULL) {
  if (is.null(treeOrder <- attr(tree, 'order')) || treeOrder != 'preorder') {
    tree <- Preorder(tree)
    if (!is.null(edgeToBreak)) {
      warning("Edge numbering modified as tree not in preorder;
               edgeToBreak and mergeEdges ignored.")
      edgeToBreak <- mergeEdges <- NULL
    }
  }
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
    if (edgeToBreak > nEdge) return(tree, TBRWarning("edgeToBreak > nEdge"))
    if (edgeToBreak < 1) return(tree, TBRWarning("edgeToBreak < 1"))
    if (edgeToBreak == 1) edgeToBreak <- which(parent == parent[1])[-1] # Use other side of root
  }
  brokenEdge <- seq_along(parent) == edgeToBreak
  brokenEdge.parentNode <- parent[edgeToBreak]
  brokenEdge.childNode  <-  child[edgeToBreak]
  
  if (!is.null(mergeEdges)) { # Quick sanity checks
    if (any(mergeEdges > nEdge)) return(TBRWarning(tree, "mergeEdges value > number of edges"))
    if (length(mergeEdges) > 2 || length(mergeEdges) == 0) 
        return(TBRWarning(tree, paste0("mergeEdges value ", paste(mergeEdges, collapse='|'),  
               " invalid; must be NULL or a vector of length 1 or 2\n  ")))
    if (length(mergeEdges) == 2 && mergeEdges[1] == mergeEdges[2]) 
      return(TBRWarning(tree, "mergeEdges values must differ"))
  }
  
  edgesCutAdrift <- DescendantEdges(edgeToBreak, parent, child, nEdge)
  edgesRemaining <- !edgesCutAdrift & !brokenEdge
  edgesOnAdriftSegment <- edgesCutAdrift | brokenEdge
  Assert(identical(edgesCutAdrift, edgesOnAdriftSegment)) # If always true, we can remove the previous line
  
  brokenEdgeParent <- child == brokenEdge.parentNode
  brokenEdgeSister <- parent == brokenEdge.parentNode & !brokenEdge
  brokenEdgeDaughters <- parent == brokenEdge.childNode
  nearBrokenEdge <- brokenEdge | brokenEdgeSister | brokenEdgeParent | brokenEdgeDaughters
  if (breakingRootEdge <- !any(brokenEdgeParent)) { 
    # Edge to break is the Root Node.
    brokenRootDaughters <- parent == child[brokenEdgeSister]
    nearBrokenEdge <- nearBrokenEdge | brokenRootDaughters
  }
  
  if (is.null(mergeEdges)) {
    candidateEdges <- which(!nearBrokenEdge & not1)
    nCandidates <- length(candidateEdges)
    if (nCandidates > 1) mergeEdges <- SampleOne(candidateEdges, len=nCandidates) else mergeEdges <- candidateEdges
  }
  if (length(mergeEdges) == 1) {
    if (edgesOnAdriftSegment[mergeEdges]) {
      adriftReconnectionEdge <- mergeEdges
      if (nearBrokenEdge[mergeEdges]) {
        samplable <- which(!edgesOnAdriftSegment & !nearBrokenEdge & not1)
      } else {
        samplable <- which(!edgesOnAdriftSegment & not1)
        if (all(edgesOnAdriftSegment == not1) && breakingRootEdge) samplable <- 1
      }
      nSamplable <- length(samplable)
      if (nSamplable == 0) return(TBRWarning(tree, "No reconnection site would modify the tree; check mergeEdge"))
      rootedReconnectionEdge <- if (nSamplable == 1) samplable else SampleOne(samplable, len=nSamplable)
      #### cat(" - Selected rooted Reconnection Edge: ", rootedReconnectionEdge, "\n")  #### DEBUGGING AID
    } else {
      rootedReconnectionEdge <- mergeEdges
      if (nearBrokenEdge[mergeEdges]) {
        samplable <- which(edgesOnAdriftSegment & !nearBrokenEdge & not1)
      } else {
        samplable <- which(edgesOnAdriftSegment & not1)
      }
      nSamplable <- length(samplable)
      if (nSamplable == 0) return(TBRWarning(tree, "No reconnection site would modify the tree; check mergeEdge"))
      adriftReconnectionEdge <- if (nSamplable == 1) samplable else SampleOne(samplable)
      #### cat(" - Selected adrift Reconnection Edge: ", adriftReconnectionEdge, "\n") #### DEBUGGING AID
    }
  } else {
    whichAdrift <- edgesOnAdriftSegment[mergeEdges]
    if (sum(whichAdrift) != 1) return(TBRWarning(tree, paste("Invalid edges selected to merge:", mergeEdges[1], mergeEdges[2])))
    adriftReconnectionEdge <- mergeEdges[whichAdrift]
    rootedReconnectionEdge <- mergeEdges[!whichAdrift]
  }
  if(nearBrokenEdge[rootedReconnectionEdge] && nearBrokenEdge[adriftReconnectionEdge]) 
    return(TBRWarning(tree, "Selected mergeEdges will not change tree topology."))
  #### edgelabels(edge = edgeToBreak, bg='orange', cex=1.8)  #### DEBUGGING AID
  #### edgelabels(edge=adriftReconnectionEdge, bg='cyan')    #### DEBUGGING AID
  #### edgelabels(edge=rootedReconnectionEdge, bg='magenta') #### DEBUGGING AID
  
  #########Assert(edgesOnAdriftSegment[adriftReconnectionEdge])
  #########Assert(!edgesOnAdriftSegment[rootedReconnectionEdge])
  
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
    #########Assert(identical(sum(repurposedDaughterEdge), sum(spareDaughterEdge), 1))
    #### which(repurposedDaughterEdge)
    #### which(spareDaughterEdge)
    child[repurposedDaughterEdge] <- child[spareDaughterEdge]
    child[spareDaughterEdge] <- parent[adriftReconnectionEdge]
    #########Assert(parent[spareDaughterEdge] == brokenEdge.childNode)
    parent[adriftReconnectionEdge] <- brokenEdge.childNode
  }
  if (!nearBrokenEdge[rootedReconnectionEdge]) {
    if (breakingRootEdge) {
      parent[brokenRootDaughters] <- brokenEdge.parentNode
      spareNode <- child[brokenEdgeSister]
      child [brokenEdgeSister] <- child[rootedReconnectionEdge]
      parent[brokenEdge | brokenEdgeSister] <- spareNode
      child[rootedReconnectionEdge] <- spareNode
    } else {
      parent[brokenEdgeSister] <- parent[brokenEdgeParent]
      parent[brokenEdgeParent] <- parent[rootedReconnectionEdge]
      parent[rootedReconnectionEdge] <- brokenEdge.parentNode
    }
  }
  
  #########Assert(identical(unique(table(parent)), 2L))
  #########Assert(identical(unique(table(child)),  1L))
  ####   matrix(c(parent, child), ncol=2)
  
  tree$edge <- RenumberTree(parent, child, nEdge)
  tree
}

#' Rooted TBR 
#' @describeIn TBR Perform \acronym{TBR} rearrangement, retaining position of root
#' @export
RootedTBR <- function(tree, edgeToBreak = NULL, mergeEdges = NULL) {
  if (is.null(treeOrder <- attr(tree, 'order')) || treeOrder != 'preorder') tree <- Preorder(tree)
  nTips <- tree$Nnode + 1
  if (nTips < 4) return (TBRWarning(tree, 'Fewer than 4 tips'))
  edge   <- tree$edge
  parent <- edge[, 1]
  child  <- edge[, 2]
  rootNode <- parent[1]
  rootEdges <- parent == rootNode
  nEdge <- length(parent)
  rightTree <- DescendantEdges(1, parent, child, nEdge)
  selectableEdges <- !rootEdges
  if (sum( rightTree) < 4) selectableEdges[ rightTree] <- FALSE
  if (sum(!rightTree) < 4) selectableEdges[!rightTree] <- FALSE
  if (!any(selectableEdges)) return(TBRWarning(tree, 'No opportunity to rearrange tree due to root position'))

  if (is.null(edgeToBreak)) {
    edgeToBreak <- SampleOne(which(selectableEdges), len=nEdge - 2L) # Pick an edge at random
  } else {
    if (edgeToBreak > nEdge) return(TBRWarning(tree, "edgeToBreak > nEdge"))
    if (edgeToBreak < 1) return(TBRWarning(tree, "edgeToBreak < 1"))
    if (rootEdges[edgeToBreak]) return(TBRWarning(tree, "RootedTBR cannot break root edge; try TBR"))
  }
  repeat {
    edgeInRight <- rightTree[edgeToBreak]
    subtreeWithRoot <- if (edgeInRight) rightTree else !rightTree
    subtreeEdges <- !rootEdges & subtreeWithRoot
    if (sum(edgesCutAdrift <- DescendantEdges(edgeToBreak, parent, child, nEdge)) > 2) break;
    if (sum(subtreeEdges, -edgesCutAdrift) > 2) break; # the edge itself, and somewheres else
    # TODO check that all expected selections are valid
    selectableEdges[edgeToBreak] <- FALSE
    ###Assert(any(selectableEdges))
    edgeToBreak <- SampleOne(which(selectableEdges), len=nEdge - 2L)
  }
  brokenEdge <- seq_along(parent) == edgeToBreak
  brokenEdge.parentNode <- parent[edgeToBreak]
  brokenEdge.childNode  <-  child[edgeToBreak]
  
  edgesRemaining <- !edgesCutAdrift & subtreeEdges
  edgesOnAdriftSegment <- edgesCutAdrift | brokenEdge
  
  if (!is.null(mergeEdges)) { # Quick sanity checks
    if (any(mergeEdges > nEdge)) return(TBRWarning(tree, "mergeEdges value > number of edges"))
    if (length(mergeEdges) > 2 || length(mergeEdges) == 0) 
        return(TBRWarning(tree, paste0("mergeEdges value ", paste(mergeEdges, collapse='|'),  
               " invalid; must be NULL or a vector of length 1 or 2\n  ")))
    if (length(mergeEdges) == 2 && mergeEdges[1] == mergeEdges[2]) 
      return(TBRWarning(tree, "mergeEdges values must differ"))
    if (!all(subtreeWithRoot[mergeEdges])) return(TBRWarning(tree, paste("mergeEdges", 
          mergeEdges[1], mergeEdges[2], "not on same side of root as edgeToBreak", edgeToBreak)))
  }  
  
  brokenEdgeParent <- child  == brokenEdge.parentNode
  brokenEdgeSister <- parent == brokenEdge.parentNode & !brokenEdge
  
  brokenEdgeDaughters <- parent == brokenEdge.childNode
  nearBrokenEdge <- brokenEdgeSister | brokenEdgeParent | brokenEdgeDaughters | brokenEdge
  ###Assert(any(brokenEdgeParent))
  
  if (is.null(mergeEdges)) {
    mergeEdges  <- which(subtreeEdges & !nearBrokenEdge)
    nCandidates <- length(mergeEdges)
    if (nCandidates > 1) mergeEdges <- SampleOne(mergeEdges, len=nCandidates)
  }
  if (length(mergeEdges) == 0) {
    return(TBRWarning(tree, paste("Breaking edge", edgeToBreak, "does not allow any new reconnections whilst preserving root position.")))
  }
  if (length(mergeEdges) == 1) {
    if (edgesOnAdriftSegment[mergeEdges]) {
      adriftReconnectionEdge <- mergeEdges
      if (nearBrokenEdge[mergeEdges]) {
        samplable <- which(subtreeEdges & !edgesOnAdriftSegment & !nearBrokenEdge)
      } else {
        samplable <- which(subtreeEdges & !edgesOnAdriftSegment)
        ###Assert(length(samplable) > 0)
      }
      nSamplable <- length(samplable)
      if (nSamplable == 0) return(TBRWarning(tree, "No reconnection site would modify the tree; check mergeEdge"))
      rootedReconnectionEdge <- if (nSamplable == 1) samplable else SampleOne(samplable, len=nSamplable)
      #### cat(" - Selected rooted Reconnection Edge: ", rootedReconnectionEdge, "\n")  #### DEBUGGING AID
    } else {
      rootedReconnectionEdge <- mergeEdges
      if (nearBrokenEdge[mergeEdges]) {
        samplable <- which(subtreeEdges & edgesOnAdriftSegment & !nearBrokenEdge)
      } else {
        samplable <- which(subtreeEdges & edgesOnAdriftSegment)
      }
      nSamplable <- length(samplable)
      if (nSamplable == 0) return(TBRWarning(tree, "No reconnection site would modify the tree; check mergeEdge"))
      adriftReconnectionEdge <- if (nSamplable == 1) samplable else SampleOne(samplable, len=nSamplable)
      #### cat(" - Selected adrift Reconnection Edge: ", adriftReconnectionEdge, "\n") #### DEBUGGING AID
    }
  } else {
    whichAdrift <- edgesOnAdriftSegment[mergeEdges]
    if (sum(whichAdrift) != 1) return(TBRWarning(tree, paste("Invalid edges selected to merge:",
            mergeEdges[1], mergeEdges[2], " - etb= ", edgeToBreak)))
    adriftReconnectionEdge <- mergeEdges[whichAdrift]
    rootedReconnectionEdge <- mergeEdges[!whichAdrift]
  }
  if(nearBrokenEdge[rootedReconnectionEdge] && nearBrokenEdge[adriftReconnectionEdge]) 
    return(TBRWarning(tree, "Selected mergeEdges will not change tree topology."))
  #### edgelabels(edge = edgeToBreak, bg='orange', cex=1.8)  #### DEBUGGING AID
  #### edgelabels(edge=adriftReconnectionEdge, bg='cyan')    #### DEBUGGING AID
  #### edgelabels(edge=rootedReconnectionEdge, bg='magenta') #### DEBUGGING AID
  
  ###Assert(edgesOnAdriftSegment[adriftReconnectionEdge])
  ###Assert(!edgesOnAdriftSegment[rootedReconnectionEdge])
  
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
    ###Assert(identical(sum(repurposedDaughterEdge), sum(spareDaughterEdge), 1))
    #### which(repurposedDaughterEdge)
    #### which(spareDaughterEdge)
    child[repurposedDaughterEdge] <- child[spareDaughterEdge]
    child[spareDaughterEdge] <- parent[adriftReconnectionEdge]
    ###Assert(parent[spareDaughterEdge] == brokenEdge.childNode)
    parent[adriftReconnectionEdge] <- child[edgeToBreak]
  }
  if (!nearBrokenEdge[rootedReconnectionEdge]) {
    parent[brokenEdgeSister] <- parent[brokenEdgeParent]
    parent[brokenEdgeParent] <- parent[rootedReconnectionEdge]
    parent[rootedReconnectionEdge] <- brokenEdge.parentNode
  }
  
  ###Assert(identical(unique(table(parent)), 2L))
  ###Assert(identical(unique(table(child)),  1L))
  ####   matrix(c(parent, child), ncol=2)
  
  retTree <- tree
  retTree$edge <- RenumberTree(parent, child, nEdge)
  retTree
}
