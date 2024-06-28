#' Collapse areas of agreement between two trees
#' 
#' ReduceTreess trees according to the tree reduction rules of
#' \insertCite{Allen2001;textual}{TreeDist}:
#' - Collapse identical pendant subtrees;
#' - Compress equivalent internal chains.
#' 
#' @template tree12Params
#' @param check Logical specifying whether to validate input. Specify
#' `FALSE` and you will encounter undefined behaviour if trees are not
#' binary `phylo` objects with identical leaf labels, rooted on leaf 1.
#' 
#' @return `ReduceTrees()` returns a list of two trees, corresponding to 
#' `tree1` and `tree2` after any identical groupings have been collapsed,
#' with tree edges listed in postorder; or `NULL` if the trees are equivalent.
#' @examples
#' tree1 <- TreeTools::BalancedTree(9)
#' tree2 <- TreeTools::PectinateTree(9)
#' 
#' # Set graphical parameters
#' oPar <- par(mai = rep(0.1, 4), mfrow = c(2, 2))
#' 
#' plot(tree1)
#' plot(tree2)
#' 
#' # Reduce trees by collapsing identical clades
#' confl <- ReduceTrees(tree1, tree2)
#' 
#' plot(confl[[1]])
#' plot(confl[[2]])
#' 
#' # Restore graphical parameters
#' par(oPar)
#' @template MRS
#' @importFrom TreeTools NTip PostorderOrder RenumberTips RootTree
#' @export
ReduceTrees <- function(tree1, tree2, check = TRUE) {
  if (check) {
    if (!inherits(tree1, 'phylo')) {
      stop("`tree1` must be a `phylo` object.")
    }
    if (!inherits(tree2, 'phylo')) {
      stop("`tree2` must be a `phylo` object.")
    }
    if (NTip(tree1) != NTip(tree2)) {
      stop("Trees must bear same leaf labels")
    }
    
    tree2 <- RenumberTips(tree2, tree1)
    tree1 <- RootTree(tree1, 1)
    tree2 <- RootTree(tree2, 1)
    if (tree1[["Nnode"]] != tree2[["Nnode"]]
        || tree1[["Nnode"]] != NTip(tree1) - 1L) {
      stop("Trees must be binary")
    }
  }
  ret <- ReduceTrees_trees(tree1$edge[PostorderOrder(tree1), ],
                      tree2$edge[PostorderOrder(tree2), ],
                      tree1[["tip.label"]])
  
  # Return:
  if (is.null(ret[[1]])) {
    NULL
  } else {
    ret
  }
}
