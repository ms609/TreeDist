# Tree-pair generators (as dev/rogue_bound/make_workload.R)

# tree1 after `nMoves` random leaf prune-and-regraft moves
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
