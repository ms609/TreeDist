#' Size of maximum agreement subtree
#' 
#' Calculate the maximum agreement subtree between two phylogenetic trees.
#' 
#' @template tree12params
#' @param tree Logical specifying whether to return the maximum agreement 
#' subtree; if `FALSE`, only its size will be returned.
#' @param rooted Logical specifying whether to treat the trees as rooted.
#' 
#' @seealso [`phangorn::mast`], a slower, all-R implementation.
#' 
#' @references 
#' \insertref{Valiente2009}{TreeDist}
#' 
#' @template MRS
#' @importFrom ape root
#' @export
MAST <- function (tree1, tree2, tree = FALSE, rooted = TRUE) {
  label1 <- tree1$tip.label
  label2 <- tree2$tip.label
  
  tree1 <- Postorder(Preorder(drop.tip(tree1, setdiff(label1, label2))))
  label1 <- tree1$tip.label
  tree2 <- Postorder(Preorder(RenumberTips(drop.tip(tree2, setdiff(label2, label1)), label1)))
  
  nTip <- length(label1)
  
  if (!rooted) {
    if (!is.rooted(tree1)) {
      tree1 <- root(tree1, outgroup = tree1$edge[nTip * 2 - 2],
                    resolve.root = TRUE)
    }
    max(vapply(seq_len(nTip - 2L) + nTip, function (node) {
      MAST(tree1, root(tree2, node=node, resolve.root = TRUE),
           tree = FALSE, rooted = TRUE)
    }, 0L))
  } else {
    if (!is.rooted(tree1) || !is.rooted(tree2)) {
      stop("Both trees must be rooted if rooted = TRUE")
    }
    cpp_mast(tree1$edge, tree2$edge, nTip)
  }
}