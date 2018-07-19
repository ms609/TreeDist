#' Frequency of splits
#' 
#' `SplitFrequency` provides a simple way to count the number of times that 
#' bipartition splits, as defined by a reference tree, occur in a forest of trees.
#' 
#' If multiple calculations are required, some time can be saved by using the
#' constituent functions (see examples)
#' 
#' 
#' @param referenceTree A tree of class phylo, or a table specifying its splits 
#'                     (as obtained through [ForestSplits])
#' @param forest a list of trees of class phylo, or a multiPhylo object; or a table
#'               enumerating the occurrences of each split in that forest 
#'                     (as obtained through [ForestSplits])
#'                     
#' @return Number of trees in `forest` that contain each split in `referenceTree`.
#'         Note that the three nodes at the root of the tree correspond to a 
#'         single split; see the example for how these might be plotted on a tree.
#' 
#' @author Martin R. Smith
#' @importFrom phangorn Descendants
#' @export
#' 
#' @examples {
#'   library(ape)
#'   set.seed(0) # Set seed so random trees are reproducable
#'   tree1 <- rtree(7)
#'   tree2 <- rtree(7)
#'   tree3 <- rtree(7)
#'   forest <- list(tree1, tree2, tree2, tree3, rtree(7))
#'   
#'   # Simple, but means counting each split in the forest twice:
#'   SplitFrequency(tree1, forest)
#'   SplitFrequency(tree2, forest)
#'   
#'   # Takes more typing, but caches the forest splits to save time at runtime:
#'   forestSplits <- ForestSplits(forest)
#'   tree1Freqs <- SplitFrequency(tree1, forestSplits)
#'   tree2Freqs <- SplitFrequency(tree2, forestSplits)
#'   
#'   plot(tree1)
#'   nodelabels(c('-', '-', # The first two nodes do not denote unique splits
#'                tree1Freqs))
#'#' }
SplitFrequency <- function(referenceTree, forest) {
  if (class(referenceTree) %in% c('list', 'phylo')) {
    referenceTree <- ForestSplits(referenceTree)
  }
  if (class(forest) %in% c('list', 'phylo', 'multiPhylo')) {
    forest <- ForestSplits(forest)
  }
  
  # Return:
  forest[names(referenceTree)]
}

#' @describeIn SplitFrequency Assign a unique integer to each split
#' @param tips Integer vector specifying the tips of the tree within the chosen split
#' @template treeParam
#' @param tipIndex Character vector of tip names, in a fixed order
#' @param powersOf2 Integer vector of same length as tipIndex, specifying a power 
#'  of 2 to be associated with each tip in turn
#' @export
SplitNumber <- function (tips, tree, tipIndex, powersOf2) {
  included <- tipIndex %in% tree$tip.label[tips]
  as.character(min(c(sum(powersOf2[included]), sum(powersOf2[!included]))))
}

#' @describeIn SplitFrequency Frequency of splits in a given forest of trees
#' @importFrom gmp as.bigz
#' @export
ForestSplits <- function (forest) {
  if (class(forest) == 'phylo') forest <- list(forest)
  tipIndex <- sort(forest[[1]]$tip.label)
  nTip <- length(tipIndex)
  powersOf2 <- as.bigz(2L ^ (seq_len(nTip) - 1L))
  splits <- table(vapply(forest, function (tr) {
    # +2: Don't consider root node (not a node) or first node (duplicated)
    vapply(Descendants(tr, nTip + 2L + seq_len(nTip - 3L), type='tips'),
           SplitNumber, character(1), tr, tipIndex, powersOf2)
  }, character(nTip - 3L)))
  # Return:
  splits[names(splits) != '0'] # 0 will occur when a tree contains polytomies
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

