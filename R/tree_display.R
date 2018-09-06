#' Consensus without taxa
#' 
#' Displays a consensus plot with selected taxa excluded.
#' 
#' A useful way to gain resolution if a few wildcard taxa obscure a consistent
#' set of relationship.
#' 
#' @param trees A list of phylogenetic trees, of class `multiPhylo` or `list`
#' @param tip A character vector specifying the names (or numbers) of tips to
#'                drop (using ape::drop.tip)
#' @param \dots Additional parameters to pass on to ape::[consensus] or [legend]
#'                
#' @return A consensus tree without the excluded taxa
#' @author Martin R. Smith
#' @importFrom ape consensus drop.tip
#' @export
ConsensusWithout <- function (trees, tip, ...) {
  if (class(trees) == 'phylo') {
    drop.tip(trees, tip=tip) 
  } else {
    consensus(lapply(trees, drop.tip, tip=tip), ...)
  }
}

#' @describeIn ConsensusWithout Adds missing taxa to a plotted consensus tree
#' @param position Where to plot the missing taxa.  See [legend] for options.
#' @importFrom graphics legend
#' @author Martin R. Smith
#' @export
MarkMissing <- function (tip, position='bottomleft', ...) {
  if (length(tip) > 0) {
    legend(position, legend=gsub('_', ' ', tip, fixed=TRUE),
         lwd=1, lty=2, bty='n', ...)
  }
}

#' Sort tree
#'
#' Sorts each node into a consistent order, so similar trees look visually similar.
#'
#' @template treeParam
#'
#' @return A tree of class phylo, with each node sorted such that the larger clade is first.
#'
#' @author Martin R. Smith
#' @export
SortTree <- function(tree) {
  edge <- tree$edge
  parent <- edge[, 1]
  child <- edge[, 2]
  tipLabels <- tree$tip.label
  tree.ntip <- length(tipLabels)
  descendants <- Descendants(tree)
  nDescendants <- vapply(descendants, length, integer(1))
  MinKid <- function (tips) min(tipLabels[tips])
  swaps <- vapply(tree.ntip + 1:tree$Nnode, function(node) {
    kids <- child[parent == node]
    descs <- nDescendants[kids]
    if (all(descs == 1L)) {
      order(tipLabels[kids])[1] == 1
    } else if (descs[1] == descs[2]) {
      order(vapply(descendants[kids], MinKid, character(1)))[1] == 1
    } else {
      descs[1] < descs[2]
    }
  }, logical(1))
  for (node in tree.ntip + rev(which(swaps))) {
    childEdges <- parent==node
    kids <- child[childEdges]
    child[childEdges][2:1] <- kids
  }
  tree$edge[, 1] <- parent
  tree$edge[, 2] <- child
  attr(tree, 'order') <- NULL
  Cladewise(Renumber(tree))
}

#' Newick Tree
#' 
#' Writes a tree in Newick format
#' 
#' @template treeParam
#' 
#' @return A character string describing `tree` in Newick format
#' 
#' @importFrom ape write.tree
#' @export
NewickTree <- function(tree) gsub('_', ' ', write.tree(tree), fixed=TRUE)
