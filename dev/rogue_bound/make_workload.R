# Generate fixed workloads of tree pairs for SPRDist(method = "Rogue").
# Random pairs: two independent RandomTree(n, root = TRUE).
# Near pairs: tree2 = tree1 after `nMoves` leaf prune-and-regraft moves.
suppressMessages(library(TreeTools, lib.loc = "/root/Rlib"))

MovedTree <- function(tree, nMoves) {
  for (i in seq_len(nMoves)) {
    leaf <- sample(tree[["tip.label"]], 1)
    pruned <- DropTip(tree, leaf)
    nEdge <- dim(pruned[["edge"]])[1]
    where <- pruned[["edge"]][sample.int(nEdge, 1), 2]
    tree <- AddTip(pruned, where = where, label = leaf)
  }
  RenumberTips(Preorder(RootTree(tree, 1)), sort(tree[["tip.label"]]))
}

RandomPair <- function(nTip) {
  t1 <- RandomTree(nTip, root = TRUE)
  t2 <- RandomTree(nTip, root = TRUE)
  list(t1, RenumberTips(t2, t1))
}

NearPair <- function(nTip, nMoves) {
  t1 <- RandomTree(nTip, root = TRUE)
  t1 <- RenumberTips(t1, sort(t1[["tip.label"]]))
  list(t1, RenumberTips(MovedTree(t1, nMoves), t1))
}

set.seed(20261014)
workloads <- list(
  rand30 = lapply(1:40, function(i) RandomPair(30)),
  rand80 = lapply(1:12, function(i) RandomPair(80)),
  rand150 = lapply(1:4, function(i) RandomPair(150)),
  near80 = lapply(1:30, function(i) NearPair(80, sample(2:6, 1))),
  near150 = lapply(1:16, function(i) NearPair(150, sample(3:8, 1)))
)
saveRDS(workloads, "workloads.rds")

# Larger set for output verification (>= 200 pairs)
set.seed(4242)
verify <- c(
  lapply(1:60, function(i) RandomPair(sample(8:40, 1))),
  lapply(1:40, function(i) RandomPair(sample(41:100, 1))),
  lapply(1:10, function(i) RandomPair(sample(101:160, 1))),
  lapply(1:60, function(i) NearPair(sample(10:60, 1), sample(1:5, 1))),
  lapply(1:40, function(i) NearPair(sample(61:150, 1), sample(2:8, 1)))
)
saveRDS(verify, "verify.rds")
cat(vapply(workloads, length, 1), length(verify), "\n")
