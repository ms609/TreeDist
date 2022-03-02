#' Collapse areas of agreement between two trees
#' 
#' Reduces trees according to the tree reduction rules of
#' \insertCite{Allen2001;textual}{TreeDist}:
#' - Collapse identical pendant subtrees;
#' - Compress equivalent internal chains.
#' 
#' @template tree12Params
#' @param check Logical specifying whether to validate input. Specify
#' `FALSE` and you will encounter undefined behaviour if trees are not
#' binary `phylo` objects with identical leaf labels, rooted on leaf 1.
#' 
#' @return `Reduce()` returns a list of two trees, corresponding to 
#' `tree1` and `tree2` after any identical groupings have been collapsed,
#' with tree edges listed in preorder; or `NULL` if the trees are equivalent.
#' @examples 
#' tree1 <- TreeTools::BalancedTree(9)
#' tree2 <- TreeTools::PectinateTree(9)
#' par(mai = rep(0.1, 4), mfrow = c(2, 2))
#' plot(tree1)
#' plot(tree2)
#' confl <- Reduce(tree1, tree2)
#' plot(confl[[1]])
#' plot(confl[[2]])
#' @template MRS
#' @importFrom TreeTools NTip PostorderOrder RenumberTips RootTree
#' @export
Reduce <- function(tree1, tree2, check = TRUE) {
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
  ret <- reduce_trees(tree1$edge[PostorderOrder(tree1), ],
                      tree2$edge[PostorderOrder(tree2), ])
  ret1 <- ret[[1]]
  if (is.null(ret1)) {
    return(NULL)
  }
  ret1Child <- ret1[, 2]
  ret1Tips <- ret1Child < 0
  tipNums <- -ret1Child[ret1Tips]
  ret1[ret1Tips, 2][sort.list(tipNums, method = "radix")] <- seq_along(tipNums)
  
  ret2 <- ret[[2]]
  ret2Child <- ret2[, 2]
  ret2Tips <- ret2Child < 0
  tipNums <- -ret2Child[ret2Tips]
  ret2[ret2Tips, 2][sort.list(tipNums, method = "radix")] <- seq_along(tipNums)
  
  keptTips <- seq_along(tree1[["tip.label"]]) %in% tipNums
  newLabs <- tree1[["tip.label"]][keptTips]
  .Retree <- function(edge) {
    structure(list(edge = edge,
                   Nnode = dim(edge)[1] / 2,
                   tip.label = newLabs),
              order = "cladewise",
              class = "phylo")
  }
  
  # Return:
  ret <- structure(list(.Retree(ret1), .Retree(ret2)), class = "multiPhylo")
  ret
  
}
