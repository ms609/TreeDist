#' Maximum Agreement Subtree size
#' 
#' Calculate the size or phylogenetic information content (Steel & Penny 2006)
#' of the maximum agreement subtree between two phylogenetic trees, i.e.
#' the largest tree that can be obtained from both `tree1` and `tree2` by
#' deleting, but not rearranging, leaves, using the algorithm of Valiente
#' (2009).
#' 
#' Implemented for trees with up to 4096 tips.  Contact maintainer if you
#' need to process larger trees.
#' 
#' @template tree12Params
#' @param rooted Logical specifying whether to treat the trees as rooted.
#' 
#' @return `MASTSize()` returns an integer specifying the number of leaves in
#' the maximum agreement subtree.
#' 
#' @examples
#' library('TreeTools') # for as.phylo, BalancedTree, PectinateTree
#' MASTSize(PectinateTree(8), BalancedTree(8))
#' MASTInfo(PectinateTree(8), BalancedTree(8))
#' 
#' MASTSize(BalancedTree(7), as.phylo(0:3, 7))
#' MASTSize(as.phylo(0:3, 7), PectinateTree(7))
#' 
#' MASTInfo(BalancedTree(7), as.phylo(0:3, 7))
#' MASTInfo(as.phylo(0:3, 7), PectinateTree(7))
#'
#' MASTSize(list(Bal = BalancedTree(7), Pec = PectinateTree(7)),
#'          as.phylo(0:3, 7))
#' MASTInfo(list(Bal = BalancedTree(7), Pec = PectinateTree(7)),
#'          as.phylo(0:3, 7))
#' 
#' CompareAll(as.phylo(0:4, 8), MASTSize)
#' CompareAll(as.phylo(0:4, 8), MASTInfo)
#' 
#' @seealso [`phangorn::mast`], a slower implementation that also lists the
#' leaves contained within the subtree.
#' 
#' @references 
#' \insertRef{Steel2006}{TreeDist}
#' 
#' \insertRef{Valiente2009}{TreeDist}
#' 
#' @template MRS
#' @family tree distances
#' @export
MASTSize <- function (tree1, tree2 = tree1, rooted = TRUE) {
  # Checks not necessary, as tip labels need not match.
  .TreeDistance(.MASTSizeSingle, tree1, tree2, rooted = rooted, checks = FALSE)
}

#' @importFrom ape drop.tip
#' @importFrom TreeTools Postorder RenumberTips TreeIsRooted RootOnNode
.MASTSizeSingle <- function (tree1, tree2, rooted = TRUE,
                             tipLabels = tree1$tip.label,
                             ...) {
  label1 <- tipLabels
  label2 <- tree2$tip.label
  
  tree1 <- drop.tip(tree1, setdiff(label1, label2))
  label1 <- tree1$tip.label
  tree2 <- RenumberTips(drop.tip(tree2, setdiff(label2, label1)), label1)
  
  nTip <- length(label1)
  
  if (!rooted) {
    if (!TreeIsRooted(tree1)) {
      tree1 <- RootOnNode(tree1, node = tree1$edge[nTip + nTip - 2L], TRUE)
    }
    postorderEdge1 <- Postorder(tree1$edge)
    tree2 <- Preorder(tree2)
    max(vapply(tree2$edge[, 2], function (node)
      .MASTSizeEdges(postorderEdge1,
                     RootOnNode(tree2, node = node, TRUE)$edge,
                     nTip = nTip), 0L))
  } else {
    if (!TreeIsRooted(tree1) || !TreeIsRooted(tree2)) {
      stop("Both trees must be rooted if rooted = TRUE")
    }
    .MASTSizeEdges(Postorder(tree1$edge), tree2$edge, nTip)
  }
}

#' @param edge1 MUST BE IN POSTORDER!
.MASTSizeEdges <- function (edge1, edge2, nTip) {
  cpp_mast(edge1 - 1L, Postorder(edge2) - 1L, nTip)
}

#' @rdname MASTSize
#' @return `MASTInfo()` returns a vector or matrix listing the phylogenetic
#' information content, in bits, of the maximum agreement subtree.
#' @importFrom TreeTools LnRooted.int LnUnrooted.int
#' @export
MASTInfo <- function (tree1, tree2 = tree1, rooted = TRUE) {
  size <- MASTSize(tree1, tree2, rooted = rooted)
  ln <- if (rooted) LnRooted.int(size) else LnUnrooted.int(size)
  ret <- ln / log(2)
  if (!is.null(attributes(size))) attributes(ret) <- attributes(size)
  # Return:
  ret
}