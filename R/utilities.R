#' Maximum splits in an _n_-leaf tree
#' 
#' `SplitsInBinaryTree()` is a convenience function to calculate the number of splits
#' in a fully-resolved (binary) tree with _n_ leaves.
#' 
#TODO delete once we can require TreeTools 1.0.1
#' This function will soon be removed from 'TreeDist' and included instead
#' in 'TreeTools' v1.0.1.
#' 
#' @param tree An object of a supported format that represents a tree or 
#' set of trees, from which the number of leaves will be calculated.
#' 
#' @return `SplitsInBinaryTree()` returns an integer vector detailing the number of 
#' unique non-trivial splits in a binary tree with _n_ leaves.
#' 
#' @examples 
#' tree <- TreeTools::BalancedTree(8)
#' SplitsInBinaryTree(tree)
#' SplitsInBinaryTree(TreeTools::as.Splits(tree))
#' SplitsInBinaryTree(8)
#' SplitsInBinaryTree(list(tree, tree))
#' @template MRS
#' @keywords internal
#' @export
SplitsInBinaryTree <- function (tree) UseMethod('SplitsInBinaryTree')

#' @rdname SplitsInBinaryTree
#' @export
#' @importFrom TreeTools NTip
SplitsInBinaryTree.phylo <- function (tree) NTip(tree) - 3L

#' @rdname SplitsInBinaryTree
#' @export
SplitsInBinaryTree.Splits <- SplitsInBinaryTree.phylo

#' @rdname SplitsInBinaryTree
#' @export
SplitsInBinaryTree.list <- function (tree) vapply(tree, SplitsInBinaryTree, integer(1L))

#' @rdname SplitsInBinaryTree
#' @export
SplitsInBinaryTree.multiPhylo <- SplitsInBinaryTree.list

#' @rdname SplitsInBinaryTree
#' @export
SplitsInBinaryTree.numeric <- function (tree) as.integer(tree) - 3L
