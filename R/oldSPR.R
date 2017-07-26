#' Subtree Pruning and Rearrangement 
#'
#' Perform one \acronym{SPR} rearrangement on a tree
#' 
#' Equivalent to phangorn's kSPR, but faster.
#' Note that rearrangements that only change the position of the root WILL be returned by 
#' \code{SPR}.  If the position of the root is irrelevant (as in Fitch parsimony, for example)
#' then this function will occasionally return a functionally equivalent topology.  
#' \code{RootIrrelevantSPR} will search tree space more efficiently in these cases.
#'
#' @template treeParam
#' @author Martin R. Smith
#' @export
oldSPR <- function(tree) { 
  tip.label <- tree$tip.label
  nTips  <- length(tip.label)
  if (nTips < 4) stop ('must be >3 tips for SPR rearrangement!')
  edge   <- tree$edge
  parent <- edge[, 1]
  child  <- edge[, 2]
  nEdge  <- length(child)
  root   <- nTips + 1
  pruning.candidates <- seq(nEdge + 1)[-root]
  repeat {
    prune.node <- SampleOne(pruning.candidates)
    moving.subnodes <- c(prune.node, which(DoDescendants(parent, child, nTips, prune.node)))
    moving.nodes <- c(prune.parent <- parent[child==prune.node], moving.subnodes)
    dont.graft.here <- c(moving.nodes, child[parent==prune.parent])
    graft.node <- c(pruning.candidates[!pruning.candidates %in% dont.graft.here])
    if (length(graft.node) > 1) graft.node <- SampleOne(graft.node)
    if (any(graft.node)) break;
    pruning.candidates <- pruning.candidates[-match(prune.node, pruning.candidates)]
    if (!any(pruning.candidates)) stop('No place to graft pruned tree')
  } 
  
  graft.edge   <- match(graft.node, child)
  graft.parent <- parent[graft.edge]
  graft.child  <-  child[graft.edge]
  prune.edge   <- match(prune.node, child)
  parent.duplicate <- parent
  parent.duplicate[prune.edge] <- NA
  sister.edge  <- match(prune.parent, parent.duplicate)
  if (prune.parent == root) {
    new.root <- child[parent==root]
    new.root <- new.root[new.root != prune.node]
    edge[sister.edge, 2] <- edge[graft.edge, 2]
    edge[graft.edge,  2] <- root
    new.root.spots <- edge==new.root
    edge[edge == root] <- new.root
    edge[new.root.spots] <- root
  } else {
    leading.edge <- match(prune.parent, child)
    edge[c(leading.edge, sister.edge, graft.edge), 2] <- edge[c(sister.edge, graft.edge, leading.edge), 2]
  }
  
  tree$edge <- OrderEdgesNumberNodes(edge[, 1], edge[, 2], nTips, nEdge)
  tree
}

#' @describeIn SPR Root position irrelevant SPR
#' @export
RootIrrelevantSPR <- function(tree) { 
  tip.label <- tree$tip.label
  nTips  <- length(tip.label)
  if (nTips < 4) stop ('must be >3 tips for SPR rearrangement!')
  edge   <- tree$edge
  parent <- edge[, 1]
  child  <- edge[, 2]
  nEdge  <- length(child)
  root   <- nTips + 1
  pruning.candidates <- seq(nEdge + 1)[-root]
  repeat {
    prune.node <- SampleOne(pruning.candidates)
    moving.subnodes <- c(prune.node, which(DoDescendants(parent, child, nTips, prune.node)))
    moving.nodes <- c(prune.parent <- parent[child==prune.node], moving.subnodes)
    if (prune.parent == root) {
      rootChildren <- child[parent==root]
      dont.graft.here <- c(moving.nodes, rootChildren, child[parent %in% rootChildren])
    } else {
      dont.graft.here <- c(moving.nodes, child[parent==prune.parent])
    }
    graft.node <- c(pruning.candidates[!pruning.candidates %in% dont.graft.here])
    if (length(graft.node) > 1) graft.node <- SampleOne(graft.node)
    if (any(graft.node)) break;
    pruning.candidates <- pruning.candidates[-match(prune.node, pruning.candidates)]
    if (!any(pruning.candidates)) stop('No place to graft pruned tree')
  } 
  
  graft.edge   <- match(graft.node, child)
  graft.parent <- parent[graft.edge]
  graft.child  <-  child[graft.edge]
  prune.edge   <- match(prune.node, child)
  parent.duplicate <- parent
  parent.duplicate[prune.edge] <- NA
  sister.edge  <- match(prune.parent, parent.duplicate)
  if (prune.parent == root) {
    new.root <- child[parent==root]
    new.root <- new.root[new.root != prune.node]
    edge[sister.edge, 2] <- edge[graft.edge, 2]
    edge[graft.edge,  2] <- root
    new.root.spots <- edge==new.root
    edge[edge == root] <- new.root
    edge[new.root.spots] <- root
  } else {
    leading.edge <- match(prune.parent, child)
    edge[c(leading.edge, sister.edge, graft.edge), 2] <- edge[c(sister.edge, graft.edge, leading.edge), 2]
  }
  
  tree$edge <- OrderEdgesNumberNodes(edge[, 1], edge[, 2], nTips, nEdge)
  tree
}

# interestingTree <- read.tree(text='(((a, (b, c)), (d, e)), (f, g));')
# interestingPruneNode <- 9
# interestingGraftNode <- 6

#' Rooted SPR rearrangement
#'
#' @importFrom ape is.rooted 
#' @importFrom stats runif 
#' @describeIn SPR Perform \acronym{SPR} operation, retaining position of root
#' @export
oldRootedSPR <- function(tree) {
  if (!is.rooted(tree)) warning("Tree root is not resolved.  Try:  tree <- SetOutgroup(tree, outgroup).")
  tip.label <- tree$tip.label
  nTips <- length(tip.label)
  if (nTips < 4) stop("Must be >3 tips for SPR rearrangement.")
  edge <- tree$edge
  parent <- edge[, 1]
  child  <- edge[, 2]
  nEdge  <- length(child)
  root   <- nTips + 1
  
  root.children <- child[parent==root]
  left.nodes <- DoDescendants(parent, child, nTips, root.children[1L])
  right.nodes <- !left.nodes
  left.nodes[root.children] <- right.nodes[root.children] <- right.nodes[root] <- FALSE
  size <- c(sum(left.nodes), sum(right.nodes))
  moves <- (size-2L) * (size-1L) / 2
  moves[size < 3] <- 0
  if (!max(moves)) return (tree)
  
  choose.right <- runif(1, min=0, max=sum(moves)) > moves[1]
  pruning.candidates <- if (choose.right) which(right.nodes) else which(left.nodes)
  subtree.base <- child[parent==root.children[choose.right + 1L]]
  subtree.basal.tip <- subtree.base < root
  if (any(subtree.basal.tip)) pruning.candidates <- pruning.candidates[-match(subtree.base[!subtree.basal.tip], pruning.candidates)]
  
  prune.node <- SampleOne(pruning.candidates)
  moving.subnodes <- c(prune.node, which(DoDescendants(parent, child, nTips, prune.node)))
  moving.nodes <- c(prune.parent <- parent[child==prune.node], moving.subnodes)
  dont.graft.here <- c(moving.nodes, child[parent==prune.parent])
  graft.candidates <- c(root.children[choose.right + 1L], pruning.candidates)
  graft.candidates <- graft.candidates[!graft.candidates %in% dont.graft.here]

  graft.child  <- SampleOne(graft.candidates)
  graft.edge   <- match(graft.child, child)
  graft.parent <- parent[graft.edge]
  
  leading.edge <- match(prune.parent, child)
  prune.edge <- match(prune.node, child)
  parent.duplicate <- parent
  parent.duplicate[prune.edge] <- NA
  sister.edge <- match(prune.parent, parent.duplicate)
  edge[c(leading.edge, sister.edge, graft.edge), 2] <- edge[c(sister.edge, graft.edge, leading.edge), 2]
  
  tree$edge <- OrderEdgesNumberNodes(edge[, 1], edge[, 2], nTips, nEdge)
  tree
}

