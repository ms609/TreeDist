#' Tree2Splits
#' 
#' Converts a phylogenetic tree to an array of bipartition splits.
#' 
#' @param tr A tree of class \code{\link[ape:read.tree]{phylo}}, with tips 
#' bearing integer labels (i.e. `tr$tip.label == 1:N`).
#' @return Returns a two-dimensional array.  Columns correspond to unique
#'  bipartitions, named with the number of a node that denotes the partition.
#'  Rows correspond to tips `1:N`.
#'
#' @author Martin R. Smith
#' 
#' @examples Tree2Splits(ape::rtree(6, tip.label=1:6, br=NULL))
#'
#' @importFrom ape reorder.phylo
#' @export
Tree2Splits <- function (tr) {
  tr <- reorder.phylo(tr, 'postorder')
  tip_label <- tr$tip.label
  n_tip <- as.integer(length(tip_label))
  root <- length(tip_label) + 1L
  bipartitions <- phangorn_bipCPP(tr$edge, n_tip)
  ret <- vapply(bipartitions[-seq_len(root)], 
                function (x) seq_len(n_tip) %in% x, 
                logical(n_tip))[seq_len(n_tip), , drop=FALSE]
  rownames(ret) <- tip_label
  colnames(ret) <- seq_len(ncol(ret)) + root
  
  ret <- UniqueSplits(ret)
  # Return:
  DropSingleSplits(ret)
}
#' @rdname Tree2Splits
#' @export
#' @keywords internal
Tree2Bipartitions <- Tree2Splits


#' Unique Splits
#' 
#' Removes equivalent duplicates from a matrix of bipartitions.
#' 
#' @param splits A logical matrix containing one named row corresponding to each
#' terminal leaf of a tree, and each column corresponds to a bipartition split;
#' each split divides terminals into two bipartitions; members of one
#' are marked `TRUE` and members of the other are marked `FALSE`.
#' @param preserveParity Logical specifying whether to preserve the `TRUE` and
#'  `FALSE` status within each split (which takes marginally longer).  If 
#'  `FALSE`, each split will be defined such that taxa in the same partition
#'  as the first element are marked `FALSE`, and other taxa marked `TRUE`.
#'  
#' @return The splits element, with all duplicate splits removed.
#' 
#' @examples 
#'   set.seed(1)
#'   splits <- Tree2Splits(ape::rtree(6, br=NULL))
#'   UniqueSplits(splits, preserveParity=TRUE)
#' 
#' @author Martin R. Smith
#' @export 
UniqueSplits <- function (splits, preserveParity = FALSE) {
  originalParity <- splits[1, ]
  
  # Set all splits to the same parity
  splits[, originalParity] <- !splits[, originalParity]
  
  # Identify duplicates
  duplicates <- duplicated(t(splits))
  
  # Remove duplicates
  ret <- splits[, !duplicates, drop=FALSE]
  
  # Return: 
  if (preserveParity) {
    ret[, originalParity[!duplicates]] <- !ret[, originalParity[!duplicates]]
  }
  
  ret
}

#' Drop Single Splits
#' 
#' Removes splits that pertain only to a single taxon from a splits object.
#' 
#' Bipartition splits are divisions, implied by each edge or node of an unrooted
#' tree topology, that divide the taxa into two groups (one of which is a clade).
#' 
#' By default, a list of splits will include those that separate a single taxon
#' (a leaf) from all others.  Such splits are, by definition, present in all 
#' trees that contain that taxon; they are not of interest when comparing trees.
#' This function removes such splits from a list of bipartitions.
#' 
#' @param split A matrix in which each column corresponds to a bipartition split
#' 
#' @return The input matrix, with any columns that separate only a single pendant
#'  tip removed.
#'         
#' @author Martin R. Smith
#' 
#' @export 
DropSingleSplits <- function (split) {
  split[, colSums(split) > 1 & colSums(!split) > 1, drop=FALSE]
}
