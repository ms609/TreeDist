#' Edge motion calcuation
#' 
#' `EdgeMotion()` measures the distance moved by each edge in a rooted
#' tree, where relative splitting times matter.
#' An edge is defined by its descendant node.
#' 
#' @template tree12Params
#' 
#' @examples
#' library("TreeTools")
#' tree1 <- as.phylo(10, 10)
#' tree2 <- as.phylo(20, 10)
#' tree1[["edge.length"]] <- rep(1, 18)
#' tree2[["edge.length"]] <- rep(1, 18)
#' 
#' tree1 <- treeB
#' tree2 <- treeA
#' tree2 <- RenumberTips(tree2, tree1)
#' tipMotion <- colSums(TipMotion(tree1, tree2))
#' maxMotion <- max(tipMotion)
#' 
#' # Configure plotting area
#' oPar <- par(mfrow = c(1, 2), mar = rep(0.1, 4))
#' 
#' # Plot trees
#' plot(tree1)
#' tiplabels(bg = hcl.colors(maxMotion + 1)[tipMotion + 1])
#' plot(tree2)
#' tiplabels(bg = hcl.colors(maxMotion + 1)[tipMotion + 1])
#' # Add legend
#' PlotTools::SpectrumLegend(
#'   "topleft",
#'   palette = hcl.colors(maxMotion + 1),
#'   title = "Tip Motion",
#'   legend = round(seq(maxMotion, 0, len = 4), 1),
#'   bty = "n")
#' 
#' # Restore plotting parameters
#' par(oPar)
#' @template MRS
#' @importFrom TreeTools RenumberTips
#' @export
TipMotion <- function(tree1, tree2) {
  heritage1 <- Heritage(tree1)
  heritage2 <- Heritage(RenumberTips(tree2, tree1))
  diff <- heritage1 - heritage2
  # Return
  sqrt(diff * diff)
}

#' @importFrom TreeTools DescendantTipsNTip PostorderOrder
#' @returns `Heritage()` returns a symmetrical matrix in which each entry
#' describes the amount of evolutionary history shared between two leaves within
#' `tree`.
Heritage <- function(tree) {
  nTip <- NTip(tree)
  nNode <- tree[["Nnode"]] + nTip
  edge <- tree[["edge"]]
  parent <- edge[, 1]
  child <- edge[, 2]
  desc <- DescendantTips(parent, child, nEdge = nEdge, nTip = nTip)
  
  length <- tree[["edge.length"]]
  if (is.null(length)) {
    length <- rep(1, length(parent))
  }
  heritage <- matrix(0, nNode, nTip)
  for (edg in rev(PostorderOrder(tree))) {
    chi <- child[edg]
    heritage[chi, ] <- heritage[parent[edg], ]
    heritage[chi, desc[edg, ]] <- heritage[chi, desc[edg, ]] + length[edg]
  }
  
  # Return:
  heritage[seq_len(nTip), ]
}
