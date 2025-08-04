library("TreeTools")
geneTree1 <- ape::read.tree(text = "((E, D), (A3, ((C, A2), (B, A1))D2 )D1);")
plot(geneTree1, show.node.label = TRUE)
nodelabels()
edgelabels()
mtext("Gene tree 1", 1)


tree <- geneTree1
#' @importFrom TreeTools NTip
DuplicationTrees <- function(tree, duplications) {
  
  if (missing(duplications)) {
    duplications <- NTip(tree) + which(tree[["node.label"]] != "")
  }
  if (length(duplications) == 0) {
    tree
  } else {
    
    edge <- tree[["edge"]]
    if (dim(edge)[[1]] != NTip(tree) * 2 - 2) {
      stop("`tree` must be bifurcating; try `MakeTreeBinary(tree)`")
    }
    parent <- edge[, 1]
    child <- edge[, 2]
    children <- child[parent == duplications[[1]]]
    do.call(c, lapply(children, function(ch) DuplicationTrees(DropTip(tree, ch))))
  }
}
DuplicationTrees(geneTree1)
