#' Rearrange edges of a phylogenetic tree
#' 
#' Rearranges a matrix that corresponds to the edges of a phylogenetic tree,
#' returning the score of the new tree.  Will generally be called from
#' within a tree search function.
#' 
#' @details \code{RearrangeTree} performs one tree rearrangement of a
#'  specified type, and returns the score of the tree (with the given dataset).
#'  It also reports the number of times that this score was hit in the 
#'  current function call.
#' 
#' @template treeParent
#' @template treeChild
#' @param dataset Third argument to pass to \code{TreeScorer}.
#' @template treeScorerParam
#' @param scoreToBeat Double giving score of input tree.
#' @param hits Integer giving number of times the input tree has already been hit.
#' @template EdgeSwapperParam
## @param  minScore trees longer than \code{minScore}, probably the score of the best previously known tree,
##     will be discarded;
## @param returnSingle returns all trees if \kbd{FALSE} or a randomly selected tree if \kbd{TRUE}.
#' @param iter iteration number of calling function, for reporting to user only.
## @template clusterParam
#' @template verbosityParam
#' @template treeScorerDots
#'
#' @author Martin R. Smith
#'
#' @template returnEdgeList
#' 
#' @examples
#' data('Lobo')
#' random.tree <- RandomTree(Lobo.phy)
#' edge <- random.tree$edge
#' parent <- edge[, 1]
#' child <- edge[, 2]
#' dataset <- PhyDat2Morphy(Lobo.phy)
#' RearrangeEdges(parent, child, dataset, EdgeSwapper=RootedNNISwap)
#' 
#' @export
RearrangeEdges <- function (parent, child, dataset, TreeScorer = MorphyLength,
                            EdgeSwapper, scoreToBeat=TreeScorer(parent, child, dataset, ...),
                            iter='?', hits=0L, verbosity=0L, ...) {
  eps <- 1e-08
  rearrangedEdges <- EdgeSwapper(parent, child)
  if (class(rearrangedEdges[[1]]) == 'list') {
    # Then we've been sent a list of possible trees
    candidateScores <- vapply(rearrangedEdges, function (edges) TreeScorer(edges[[1]], edges[[2]], dataset, ...), double(1))
    candidateScore <- min(candidateScores)
    best <- candidateScores == candidateScore
    nBest <- sum(best)
    if (candidateScore > (scoreToBeat + eps)) {
      if (verbosity > 3L) cat("\n    . Iteration", iter, '- Rearranged tree score', candidateScore, 
                              "> target", scoreToBeat)
    } else if (candidateScore + eps > scoreToBeat) { # i.e. scores are equal
      hits <- hits + nBest
      if (verbosity > 2L) cat("\n    - Iteration", iter, "- Best score", scoreToBeat, 
                              "found again", nBest, "times; now found", hits, "times.")
    } else {
      hits <- nBest
      if (verbosity > 1L) cat("\n    * Iteration", iter, "- New best score", candidateScore, 
                              "found on", hits, "trees.")
    }
    rearrangedEdges <- rearrangedEdges[[SampleOne(which(best), nBest)]]
  } else {
    candidateScore <- TreeScorer(rearrangedEdges[[1]], rearrangedEdges[[2]], dataset, ...)
    if (candidateScore > (scoreToBeat + eps)) {
      if (verbosity > 3L) cat("\n    . Iteration", iter, '- Rearranged tree score', candidateScore, "> target", scoreToBeat)
    } else if (candidateScore + eps > scoreToBeat) { # i.e. scores are equal
      hits <- hits + 1L
      if (verbosity > 2L) cat("\n    - Iteration", iter, "- Best score", scoreToBeat, "hit", hits, "times.")
    } else {
      hits <- 1L
      if (verbosity > 1L) cat("\n    * Iteration", iter, "- New best score", candidateScore, "found on", hits, "trees.")
    }
  }
  rearrangedEdges[3:4] <- c(candidateScore, hits)
  # Return:
  rearrangedEdges
}

#' Root Tree on specified tips
#' 
#' Roots a tree on the smallest clade containing the specified tips.
#' 
#' @template treeParam
#' @template outgroupTipsParam
#'                     
#' @return A tree of class phylo, rooted on the smallest clade that contains the specified tips
#' 
#' @author Martin R. Smith
#' @importFrom phangorn Ancestors Descendants
#' @importFrom ape root
#' @export
RootTree <- function (tree, outgroupTips) {
  tipLabels <- tree$tip.label
  if (!all(outgroupTips %in% tipLabels)) {
    stop("Outgroup tips", paste(outgroupTips, collapse=', '), 
         "not found in tree's tip labels.")
  }
  if (length(outgroupTips) == 1) {
    outgroup <- outgroupTips
  } else {
    tipNos <- which(tipLabels %in% outgroupTips)
    ancestry <- unlist(Ancestors(tree, tipNos))
    ancestryTable <- table(ancestry)
    lineage <- as.integer(names(ancestryTable))
    lca <- max(lineage[ancestryTable == length(outgroupTips)])
    rootNode <- length(tipLabels) + 1L
    if (lca == rootNode) {
      lca <- lineage[lineage - c(lineage[-1], 0) != -1][1] + 1L
    }
    outgroup <- Descendants(tree, lca)[[1]]
  }
  
  Renumber(root(tree, outgroup, resolve.root = TRUE))
}

#' Collapse nodes on a phylogenetic tree
#' 
#' Collapses specified nodes or edges on a phylogenetic tree, resulting in
#' polytomies.
#' 
#' @template treeParam
#' @param nodes,edges Integer vector specifying the nodes or edges in the tree
#'  to be dropped. 
#' (Use \code{\link[ape]{nodelabels}} or \code{\link[ape]{edgelabels}} 
#' to view numbers on a plotted tree.)
#' 
#' @return `tree`, with the specified nodes or edges collapsed.  
#' The length of each dropped edge will (naively) be added to each descendant edge.
#' 
#' @examples 
#'   library(ape)
#'   set.seed(1)
#'   
#'   tree <- rtree(7)
#'   par(mfrow=c(2, 1), mar=rep(0.5, 4))
#'   plot(tree)
#'   nodelabels()
#'   edgelabels(round(tree$edge.length, 2), cex=0.6, frame='n', adj=c(1, -1))
#'   
#'   newTree <- CollapseNode(tree, c(12, 13))
#'   plot(newTree)
#'   nodelabels()
#'   edgelabels(round(newTree$edge.length, 2), cex=0.6, frame='n', adj=c(1, -1))
#' 
#' @author  Martin R. Smith
#' @export
CollapseNode <- function (tree, nodes) {
  if (length(nodes) == 0) return (tree)
  
  edge <- tree$edge
  lengths <- tree$edge.length
  hasLengths <- !is.null(lengths)
  parent <- edge[, 1]
  child <- edge[, 2]
  root <- min(parent)
  nTip <- root - 1L
  maxNode <- max(parent)
  edgeBelow <- order(child)
  edgeBelow <- c(edgeBelow[1:(root-1L)], NA, edgeBelow[-(1:root-1L)])
  nodes <- unique(nodes)
  
  if (class(tree) != 'phylo') stop ("tree must be an object of class phylo")
  if (!all(nodes %in% (root + 1L):maxNode)) stop("nodes must be integers between ",
                                                 root + 1L, " and ", maxNode)
  
  keptEdges <- -edgeBelow[nodes]

  for (node in rev(sort(nodes))) {
    newParent <- parent[edgeBelow[node]]
    if (hasLengths) lengths[parent == node] <- lengths[parent == node] + lengths[child == node]
    parent[parent == node] <- newParent
  }
  
  newNumber <- c(seq_len(nTip), nTip + cumsum(root:maxNode %in% parent))
  
  tree$edge <-cbind(newNumber[parent[keptEdges]], newNumber[child[keptEdges]])
  tree$edge.length <- lengths[keptEdges]
  tree$Nnode <- tree$Nnode - length(nodes)
  
  # TODO renumber nodes sequentially
  tree
}

#' @rdname CollapseNode
#' @export
CollapseEdge <- function (tree, edges) {
  CollapseNode(tree, tree$edge[edges, 2])
}

#' Check that all nodes in a tree are bifurcating.
#' 
#' @template treeParent
#' 
#' @return Returns `NULL`, but will `stop` with an error message if a tree
#' does not appear to be bifurcating.
#' 
#' @author Martin R. Smith
#' @keywords internal
#' @export
StopUnlessBifurcating <- function (parent) {
  if (!all(table(parent) == 2L)) stop ("Tree must be bifurcating; try collapse.singles or multi2di.")
}
