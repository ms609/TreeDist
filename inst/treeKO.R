library("TreeTools")
geneTree1 <- ape::read.tree(text = "(root, ((E, D), (A3, ((C, A2), (B, A1))D2 )D1));")
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

DenumberTips <- function(trees) {
  lapply(trees, function(tr) {
    tr[["tip.label"]] <- gsub("\\d+$", "", tr[["tip.label"]], perl = TRUE)
    tr
  })
}

dupTrees <- DenumberTips(DuplicationTrees(geneTree1))

par(mfrow = c(1, 3), mar = rep(0, 4))
for (tr in dupTrees) plot(tr)

# All trees have distance zero - not very interesting for next step
distMat <- ClusteringInfoDistance(dupTrees)

# Some trees that differ
dupTrees <- c(as.phylo(0:2, 8), as.phylo(2:4, 6))
distMat <- as.matrix(ClusteringInfoDistance(dupTrees, normalize = TRUE))
leavesRemaining <- NTip(dupTrees)
leavesPruned <- 8 - leavesRemaining
r <- outer(leavesRemaining, leavesRemaining, `+`)
p <- outer(leavesPruned, leavesPruned, `+`)

d <- ((distMat * r) + p) / (r + p)
pairing <- LAPJV(`diag<-`(d, 1))[["matching"]]
