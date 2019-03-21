#' Frequency of splits
#' 
#' `SplitFrequency` provides a simple way to count the number of times that 
#' bipartition splits, as defined by a reference tree, occur in a forest of trees.
#' 
#' If multiple calculations are required, some time can be saved by using the
#' constituent functions (see examples)
#' 
#' 
#' @param reference A tree of class phylo, or a character vector specifying its splits 
#'                     (as obtained through [Tree2Splits])
#' @param forest a list of trees of class phylo, or a multiPhylo object; or a
#' list of their constituent splits (as obtained through [Tree2Splits])
#'                     
#' @return Number of trees in `forest` that contain each split in `reference`. 
#'         if `reference` is a tree of class phylo, then the sequence will
#'         correspond to the order of nodes (use `ape::nodelabels to view).
#'         Note that the three nodes at the root of the tree correspond to a 
#'         single split; see the example for how these might be plotted on a tree.
#' 
#' 
#' @examples {
#'   library(ape) # for functions rtree & nodelabels
#'   set.seed(0) # Set seed so random trees are reproducable
#'   tree1 <- rtree(7)
#'   tree2 <- rtree(7)
#'   tree3 <- rtree(7)
#'   forest <- list(tree1, tree2, tree2, tree3, rtree(7))
#'   
#'   # Simple, but means counting each split in the forest twice:
#'   tree1Freqs <- SplitFrequency(tree1, forest)
#'   SplitFrequency(tree2, forest)
#'   
#'   plot(tree1)
#'   nodelabels(tree1Freqs, node=as.integer(names(tree1Freqs)))
#' }
#' 
#' @author Martin R. Smith
#' @export
SplitFrequency <- function(reference, forest) {
  if (class(reference) %in% c('list', 'phylo')) {
    referenceSplits <- Tree2Splits(reference)
  } else {
    referenceSplits <- reference
  }
  if (class(forest) == 'phylo') {
    forestSplits <- Tree2Splits(forest) 
  } else if (class(forest) %in% c('list', 'multiPhylo')) {
    forestSplits <- lapply(forest, Tree2Splits)
  } else {
    forestSplits <- forest
  }
  
  ret <- rowSums(vapply(forestSplits, function (cf) SplitsRepeated(referenceSplits, cf),
         logical(ncol(referenceSplits))))
  names(ret) <- colnames(referenceSplits)
  # Return:
  ret
}

#' Are Splits Repeated?
#' 
#' Determines whether any of the splits described in `original` are present in
#' `matches`.
#' 
#' @param original,matches Logical matrix describing the bipartition splits
#' defined by a tree (see [Tree2Splits] for format).  The first row
#' of `original` must be `FALSE` for all splits.
#' 
#' @keywords internal
#' @export
#' @author Martin R. Smith
#' @return `SplitsRepeated` returns a logical vector of length `ncol(original)`,
#' specifying whether each split listed in `original` also appears in `matches`.
SplitsRepeated <- function (original, matches) {
  comparison <- matches[rownames(original), , drop=FALSE]
  comparison[, comparison[1, ]] <- !comparison[, comparison[1, ]]
  duplicated(t(cbind(comparison, original)))[-seq_len(ncol(matches))]
}

#' @describeIn SplitFrequency Assign a unique integer to each split
#' @param tips Integer vector specifying the tips of the tree within the chosen split
#' @template treeParam
#' @param tipIndex Character vector of tip names, in a fixed order
#' @param powersOf2 Integer vector of same length as tipIndex, specifying a power 
#'  of 2 to be associated with each tip in turn
#' @export
SplitNumber <- function (tips, tree, tipIndex, powersOf2) {
  .Deprecated("SplitFrequency")
  included <- tipIndex %in% tree$tip.label[tips]
  as.character(min(c(sum(powersOf2[included]), sum(powersOf2[!included]))))
}

#' @describeIn SplitFrequency Frequency of splits in a given forest of trees
#' @export
ForestSplits <- function (forest, powersOf2) {
  .Deprecated("SplitFrequency")
  if (class(forest) == 'phylo') forest <- structure(list(forest), class='multiPhylo')
  tipIndex <- sort(forest[[1]]$tip.label)
  nTip <- length(tipIndex)
  
  # Return:
  table(vapply(forest, function (tr) {
    # +2: Don't consider root node (not a node) or first node (duplicated)
    vapply(Descendants(tr, nTip + 2L + seq_len(nTip - 3L), type='tips'),
           SplitNumber, character(1), tr, tipIndex, powersOf2)
  }, character(nTip - 3L)))
}

#' @describeIn SplitFrequency Deprecated. Listed the splits in a given tree. 
#' Use Quartet::Tree2Splits instead.
#' @export
TreeSplits <- function (tree) {
  .Deprecated("Tree2Splits")
}

#' Support colour
#' @param support A vector of doubles in the range 0-1
#' @param show1 Logical specifying whether to display values of 1. 
#'              A transparent white will be returned if `FALSE`.  
#' @return A string containing the hexadecimal code for a colour picked from a
#'         diverging scale, or `red` if a value is invalid.
#' @importFrom colorspace diverge_hcl
#' @export
SupportColour <- function (support, show1=TRUE) {
  # continuousScale <- rev(colorspace::heat_hcl(101, h=c(300, 75), c.=c(35, 95), l=c(15, 90), power=c(0.8, 1.2))) # Viridis prefered
  divergingScale <- rev(diverge_hcl(101, h=c(260, 0), c=100, l=c(50, 90), power=1.0))
  ifelse(is.na(support) | support < 0 | support > 1 | support == '', 'red',
         ifelse(support == 1 & !show1, "#ffffff00", divergingScale[(support * 100) + 1L]))
}

#' @describeIn SupportColour alternative spelling
#' @export
SupportColor <- SupportColour

