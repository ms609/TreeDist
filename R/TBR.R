#' TBR Warning
#' Print a warning and return given tree
#'
#' @param tree tree to return
#' @param error error message to report
#'
#' @return the tree specified in tree
#' @examples
#' testFunction <- function (tree) {
#'  return(TBRWarning(parent, child, 'Message text'))
#' }
#' \dontrun{testFunction(0) # will trigger warning}
#' 
#'
#' @author Martin R. Smith
#' @keywords internal
#' @export
TBRWarning <- function (parent, child, error) {
  warning ("No TBR operation performed.\n  > ", error)
  return(list(parent, child))
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
#'  \insertRef{Felsenstein2004}{TreeSearch}
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
  
  edge <- tree$edge  
  tree$edge <- ListToMatrix(TBRSwap(edge[, 1], edge[, 2], edgeToBreak=edgeToBreak, 
                                    mergeEdges=mergeEdges))
  tree
}


## TODO Do edges need to be pre-ordered before coming here?
#' @describeIn TBR faster version that takes and returns parent and child parameters
#' @template treeParent
#' @template treeChild
#' @param nEdge (optional) Number of edges.
#' @return a list containing two elements, corresponding in turn to the rearranged parent and child parameters
#' @export
TBRSwap <- function(parent, child, nEdge = length(parent), edgeToBreak=NULL, mergeEdges=NULL) {
  if (nEdge < 5) return (list(parent, child)) #TODO do we need to re-root this tree?
  
  # Pick an edge at random
  allEdges <- seq_len(nEdge - 1L) + 1L # Only include one root edge
  not1 <- !logical(nEdge)
  not1[1] <- FALSE
  if (is.null(edgeToBreak)) {
    edgeToBreak <- SampleOne(allEdges, len=nEdge - 1L)
  } else {
    if (edgeToBreak > nEdge) return(TBRWarning(parent, child, "edgeToBreak > nEdge"))
    if (edgeToBreak < 1) return(TBRWarning(parent, child, "edgeToBreak < 1"))
    if (edgeToBreak == 1) edgeToBreak <- which(parent == parent[1])[-1] # Use other side of root
  }
  brokenEdge <- seq_along(parent) == edgeToBreak
  brokenEdge.parentNode <- parent[edgeToBreak]
  brokenEdge.childNode  <-  child[edgeToBreak]
  
  if (!is.null(mergeEdges)) { # Quick sanity checks
    if (any(mergeEdges > nEdge)) return(TBRWarning(parent, child, "mergeEdges value > number of edges"))
    if (length(mergeEdges) > 2 || length(mergeEdges) == 0) 
      return(TBRWarning(parent, child, paste0("mergeEdges value ", paste(mergeEdges, collapse='|'),  
                                              " invalid; must be NULL or a vector of length 1 or 2\n  ")))
    if (length(mergeEdges) == 2 && mergeEdges[1] == mergeEdges[2]) 
      return(TBRWarning(parent, child, "mergeEdges values must differ"))
  }
  
  edgesCutAdrift <- DescendantEdges(edgeToBreak, parent, child, nEdge)
  edgesRemaining <- !edgesCutAdrift & !brokenEdge
  
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
    if (edgesCutAdrift[mergeEdges]) {
      adriftReconnectionEdge <- mergeEdges
      if (nearBrokenEdge[mergeEdges]) {
        samplable <- which(!edgesCutAdrift & !nearBrokenEdge & not1)
      } else {
        samplable <- which(!edgesCutAdrift & not1)
        if (all(edgesCutAdrift == not1) && breakingRootEdge) samplable <- 1
      }
      nSamplable <- length(samplable)
      if (nSamplable == 0) return(TBRWarning(parent, child, "No reconnection site would modify the tree; check mergeEdge"))
      rootedReconnectionEdge <- if (nSamplable == 1) samplable else SampleOne(samplable, len=nSamplable)
      #### cat(" - Selected rooted Reconnection Edge: ", rootedReconnectionEdge, "\n")  #### DEBUGGING AID
    } else {
      rootedReconnectionEdge <- mergeEdges
      if (nearBrokenEdge[mergeEdges]) {
        samplable <- which(edgesCutAdrift & !nearBrokenEdge & not1)
      } else {
        samplable <- which(edgesCutAdrift & not1)
      }
      nSamplable <- length(samplable)
      if (nSamplable == 0) return(TBRWarning(parent, child, "No reconnection site would modify the tree; check mergeEdge"))
      adriftReconnectionEdge <- if (nSamplable == 1) samplable else SampleOne(samplable)
      #### cat(" - Selected adrift Reconnection Edge: ", adriftReconnectionEdge, "\n") #### DEBUGGING AID
    }
  } else {
    whichAdrift <- edgesCutAdrift[mergeEdges]
    if (sum(whichAdrift) != 1) return(TBRWarning(parent, child, paste("Invalid edges selected to merge:", mergeEdges[1], mergeEdges[2])))
    adriftReconnectionEdge <- mergeEdges[whichAdrift]
    rootedReconnectionEdge <- mergeEdges[!whichAdrift]
  }
  if(nearBrokenEdge[rootedReconnectionEdge] && nearBrokenEdge[adriftReconnectionEdge]) 
    return(TBRWarning(parent, child, "Selected mergeEdges will not change tree topology."))
  #### edgelabels(edge = edgeToBreak, bg='orange', cex=1.8)  #### DEBUGGING AID
  #### edgelabels(edge=adriftReconnectionEdge, bg='cyan')    #### DEBUGGING AID
  #### edgelabels(edge=rootedReconnectionEdge, bg='magenta') #### DEBUGGING AID
  
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
  return (RenumberEdges(parent, child, nEdge))
}

#' @describeIn TBR Possible TBR moves
#' @param avoid Integer vector specifying which edges should not be broken
#' @param retainRoot logical specifying whether taxa may be swapped across the root
#' @return a matrix with two columns, each row listing an edge that can be broken
#'         and an edge into which it can be merged
#' @export
TBRMoves <- function(parent, child, nEdge = length(parent), avoid=NULL, retainRoot=FALSE) {
  if (nEdge < 5) stop("No TBR rearrangements possible on a tree with < 5 edges")
  
  # Pick an edge at random
  allEdges <- seq_len(nEdge - 1L) + 1L # Only include one root edge
  logicals <- diag(nEdge) == 1
  isBreakable <- !logicals[1, ]
  descendants <- AllDescendantEdges(parent, child, nEdge)
  # not1 <- isBreakable
  
  if (length(avoid) > 0 && any(avoid > nEdge || avoid < 1)) {
    warning("Invalid edges in `avoid` parameter: edges are numbered from 1 to ", nEdge)
  }
  
  if (retainRoot) {
    leftRootEdge <- match(TRUE, parent[-1] == parent[1])+ 1L
    avoid <- c(avoid, leftRootEdge)
    rightEdges <- descendants[1, ]
    leftEdges <- !rightEdges
    if (sum(rightEdges) < 4) avoid <- c(avoid, which(rightEdges))
    if (sum(leftEdges) < 4) avoid <- c(avoid, which(leftEdges))
  }

  isBreakable[avoid] <- FALSE
  if (!any(isBreakable)) return (NULL) # no rearrangements possible
  breakable <- which(isBreakable)
  
  mergeable <- lapply(breakable, function (edgeToBreak) {
    brokenEdge <- logicals[edgeToBreak, ]
    brokenEdge.parentNode <- parent[edgeToBreak]
    brokenEdge.childNode  <-  child[edgeToBreak]
  
    edgesCutAdrift <- descendants[edgeToBreak, ]
    edgesRemaining <- !edgesCutAdrift & !brokenEdge
    
    brokenEdgeParent <- child == brokenEdge.parentNode
    brokenEdgeSister <- parent == brokenEdge.parentNode & !brokenEdge
    brokenEdgeDaughters <- parent == brokenEdge.childNode
    nearBrokenEdge <- brokenEdge | brokenEdgeSister | brokenEdgeParent | brokenEdgeDaughters
    if (breakingRootEdge <- !any(brokenEdgeParent)) { 
      # Edge to break is the Root Node.
      brokenRootDaughters <- parent == child[brokenEdgeSister]
      nearBrokenEdge <- nearBrokenEdge | brokenRootDaughters
    }
    candidates <- !nearBrokenEdge & isBreakable
    if (retainRoot) {
      candidates <- candidates & if (rightEdges[edgeToBreak]) rightEdges else leftEdges
    }
    which(candidates)
  })
  matrix(c(rep(breakable, vapply(mergeable, length, 1L)), unlist(mergeable)), ncol=2)
}

#' @describeIn TBR All unique trees one TBR move away
#' @return a list of trees, in parent-child format
#' @export
AllTBR <- function (parent, child, nEdge = length(parent), avoid=NULL, retainRoot=FALSE) {
  moves <- TBRMoves(parent, child, nEdge=nEdge, avoid=avoid, retainRoot=retainRoot)
  newTrees <- apply(moves, 1, function (edges) {
    TBRSwap(parent, child, nEdge, edges[1], edges[2])
  })
  unique(newTrees)
}

#' Rooted TBR 
#' @describeIn TBR Perform \acronym{TBR} rearrangement, retaining position of root
#' @export
RootedTBR <- function(tree, edgeToBreak = NULL, mergeEdges = NULL) {
  if (is.null(treeOrder <- attr(tree, 'order')) || treeOrder != 'preorder') tree <- Preorder(tree)
  edge   <- tree$edge
  tree$edge <- ListToMatrix(RootedTBRSwap(edge[, 1], edge[, 2], 
                            edgeToBreak=edgeToBreak, mergeEdges=mergeEdges))
  tree
}

#' @describeIn TBR faster version that takes and returns parent and child parameters
#' @export
RootedTBRSwap <- function (parent, child, nEdge=length(parent), edgeToBreak=NULL, mergeEdges=NULL) {
  if (nEdge < 5) return (TBRWarning(parent, child, 'Fewer than 4 tips'))
  nTips <- (nEdge / 2L) + 1L
  rootNode <- parent[1]
  rootEdges <- parent == rootNode
  rightTree <- DescendantEdges(1, parent, child, nEdge)
  selectableEdges <- !rootEdges
  if (sum( rightTree) < 4) {
    selectableEdges[ rightTree] <- FALSE
  } else if (sum( rightTree) < 6) {
    rightChild <- child[1]
    rightGrandchildEdges   <- parent==rightChild
    rightGrandchildren     <- child[rightGrandchildEdges]
    rightGrandchildrenTips <- rightGrandchildren <= nTips
    selectableEdges[which(rightGrandchildEdges)[!rightGrandchildrenTips]] <- FALSE  
  }
  if (sum(!rightTree) < 4) {
    selectableEdges[!rightTree] <- FALSE
  } else if (sum(!rightTree) < 6) {
     leftChild <- child[rootEdges][2]
     leftGrandchildEdges   <- parent==leftChild
     leftGrandchildren     <- child[ leftGrandchildEdges]
     leftGrandchildrenTips <-  leftGrandchildren <= nTips
     selectableEdges[which( leftGrandchildEdges)[! leftGrandchildrenTips]] <- FALSE  
  }
  
  if (!any(selectableEdges)) return(TBRWarning(parent, child, 'No opportunity to rearrange tree due to root position'))

  if (is.null(edgeToBreak)) {
    edgeToBreak <- SampleOne(which(selectableEdges)) # Pick an edge at random
  } else {
    if (edgeToBreak > nEdge) return(TBRWarning(parent, child, "edgeToBreak > nEdge"))
    if (edgeToBreak < 1) return(TBRWarning(parent, child, "edgeToBreak < 1"))
    if (rootEdges[edgeToBreak]) return(TBRWarning(parent, child, "RootedTBR cannot break root edge; try TBR"))
    if (!selectableEdges[edgeToBreak]) return(TBRWarning(parent, child, paste("Breaking edge", edgeToBreak,
                                              "does not allow a changing reconnection")))
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
    edgeToBreak <- SampleOne(which(selectableEdges))
  }
  brokenEdge <- seq_along(parent) == edgeToBreak
  brokenEdge.parentNode <- parent[edgeToBreak]
  brokenEdge.childNode  <-  child[edgeToBreak]
  
  edgesRemaining <- !edgesCutAdrift & subtreeEdges
  edgesOnAdriftSegment <- edgesCutAdrift | brokenEdge
  
  if (!is.null(mergeEdges)) { # Quick sanity checks
    if (any(mergeEdges > nEdge)) return(TBRWarning(parent, child, "mergeEdges value > number of edges"))
    if (length(mergeEdges) > 2 || length(mergeEdges) == 0) 
        return(TBRWarning(parent, child, paste0("mergeEdges value ", paste(mergeEdges, collapse='|'),  
               " invalid; must be NULL or a vector of length 1 or 2\n  ")))
    if (length(mergeEdges) == 2 && mergeEdges[1] == mergeEdges[2]) 
      return(TBRWarning(parent, child, "mergeEdges values must differ"))
    if (!all(subtreeWithRoot[mergeEdges])) return(TBRWarning(parent, child, paste("mergeEdges", 
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
    return(TBRWarning(parent, child, paste("Breaking edge", edgeToBreak, "does not allow any new reconnections whilst preserving root position.")))
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
      if (nSamplable == 0) return(TBRWarning(parent, child, "No reconnection site would modify the tree; check mergeEdge"))
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
      if (nSamplable == 0) return(TBRWarning(parent, child, "No reconnection site would modify the tree; check mergeEdge"))
      adriftReconnectionEdge <- if (nSamplable == 1) samplable else SampleOne(samplable, len=nSamplable)
      #### cat(" - Selected adrift Reconnection Edge: ", adriftReconnectionEdge, "\n") #### DEBUGGING AID
    }
  } else {
    whichAdrift <- edgesOnAdriftSegment[mergeEdges]
    if (sum(whichAdrift) != 1) return(TBRWarning(parent, child, paste("Invalid edges selected to merge:",
            mergeEdges[1], mergeEdges[2], " - etb= ", edgeToBreak)))
    adriftReconnectionEdge <- mergeEdges[whichAdrift]
    rootedReconnectionEdge <- mergeEdges[!whichAdrift]
  }
  if(nearBrokenEdge[rootedReconnectionEdge] && nearBrokenEdge[adriftReconnectionEdge]) 
    return(TBRWarning(parent, child, "Selected mergeEdges will not change tree topology."))
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
  return (RenumberEdges(parent, child, nEdge))
}