test_that("MCITree returns one of the input trees", {
  trees <- as.phylo(24:40, 16)
  mci <- MCITree(trees)
  expect_true(any(vapply(trees, all.equal, TRUE, mci) == "TRUE"))
})

test_that("MCITree modes give valid results", {
  trees <- as.phylo(24:40, 16)
  for (mode in c("phylogenetic", "clustering", "credibility")) {
    result <- MCITree(trees, mode)
    expect_s3_class(result, "phylo")
  }
})

test_that("MCITree aliases match base modes", {
  trees <- as.phylo(24:40, 16)
  expect_equal(MCITree(trees, "spic"), MCITree(trees, "phylogenetic"))
  expect_equal(MCITree(trees, "scic"), MCITree(trees, "clustering"))
})

test_that("MCITree handles edge cases", {
  tree <- PectinateTree(8)
  expect_identical(MCITree(tree), tree)

  trees <- c(PectinateTree(8), BalancedTree(8))
  expect_identical(MCITree(trees[1]), trees[[1]])
})

test_that("MCITree rejects invalid info", {
  trees <- as.phylo(24:40, 16)
  expect_error(MCITree(trees, "invalid"), "must be")
})

test_that("MCITree credibility selects highest-support tree", {
  # 10 copies of tree A + 1 copy of tree B:
  # tree A's splits have frequency 10/11; tree B's have ~1/11.
  # The credibility tree should be one of the A copies.
  treeA <- PectinateTree(8)
  treeB <- BalancedTree(8)
  trees <- c(rep(list(treeA), 10), list(treeB))
  class(trees) <- "multiPhylo"
  mcc <- MCITree(trees, "credibility")
  expect_equal(mcc, treeA)
})

test_that("MCITree messages on ties", {
  # Identical trees should tie
  tree <- PectinateTree(8)
  trees <- rep(list(tree), 5)
  class(trees) <- "multiPhylo"
  expect_message(MCITree(trees), "tied")
})
