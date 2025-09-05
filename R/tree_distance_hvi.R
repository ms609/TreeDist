#' Hierarchical Variation of Information distance
#' 
#' Calculate the hierachicical variation of information distance
#' 
#' Explain here how the hierachical variation of information distance works
#' @export
HierachicalMutual <- function (tree1, tree2=NULL, ...) {
  treeA <- ape::write.tree(tree1)
  treeB <- ape::write.tree(tree2)
  MutualInformation <- d_n(treeA, treeB)
  MutualInformation
}

#' @export
HierachicalMutual <- HierachicalMutual