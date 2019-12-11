#' Size of maximum agreement subtree
#' 
#' Calculate the maximum agreement subtree between two phylogenetic trees, i.e.
#' the largest tree that can be obtained from either `tree1` or `tree2` solely
#' by deleting tips.
#' 
#' @template tree12params
#' @param rooted Logical specifying whether to treat the trees as rooted.
#' 
#' @return `MASTSize` returns an integer specifying the number of tips in the
#' maximum agreement subtree.
#' 
#' @examples
#' library('TreeTools')
#' MASTSize(PectinateTree(8), BalancedTree(8))
#' 
#' @seealso [`phangorn::mast`], a slower, all-R implementation that also returns
#' the tips contained within the subtree.
#' 
#' @references 
#' \insertRef{Valiente2009}{TreeDist}
#' 
#' @template MRS
#' @family tree distances
#' @importFrom ape drop.tip root
#' @importFrom TreeTools Postorder RenumberTips
#' @export
MASTSize <- function (tree1, tree2, rooted = TRUE) {
  label1 <- tree1$tip.label
  label2 <- tree2$tip.label
  
  tree1 <- Postorder(drop.tip(tree1, setdiff(label1, label2)))
  label1 <- tree1$tip.label
  tree2 <- Postorder(RenumberTips(drop.tip(tree2, 
                                           setdiff(label2, label1)), label1))
  
  nTip <- length(label1)
  
  if (!rooted) {
    if (!is.rooted(tree1)) {
      tree1 <- root(tree1, outgroup = tree1$edge[nTip * 2 - 2],
                    resolve.root = TRUE)
    }
    max(vapply(seq_len(nTip - 3L) + nTip + 2L, function (node) {
      MASTSize(tree1, root(tree2, node = node, resolve.root = TRUE), 
               rooted = TRUE)}, 0L))
  } else {
    if (!is.rooted(tree1) || !is.rooted(tree2)) {
      stop("Both trees must be rooted if rooted = TRUE")
    }
    cpp_mast(tree1$edge - 1L, tree2$edge - 1L, nTip)
  }
}