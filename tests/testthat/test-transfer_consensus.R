test_that("Identical trees return fully resolved consensus", {
  tree <- as.phylo(0, nTip = 10)
  trees <- structure(rep(list(tree), 10), class = "multiPhylo")

  tc <- TransferConsensus(trees, greedy = "first")
  expect_equal(NSplits(tc), NSplits(tree))

  tc_best <- TransferConsensus(trees, greedy = "best")
  expect_equal(NSplits(tc_best), NSplits(tree))
})

test_that("Star tree returned when no signal", {
  set.seed(6129)
  # Fully random trees with many tips — negligible split overlap
  trees <- structure(lapply(1:30, function(i) rtree(50)),
                     class = "multiPhylo")
  tc <- TransferConsensus(trees, greedy = "first")

  # Should be very unresolved (0 or near-0 splits)
  expect_lte(NSplits(tc), 5)
})

test_that("Transfer consensus is at least as resolved as majority-rule for structured trees",
{
  # Trees from as.phylo with moderate overlap
  trees <- as.phylo(0:29, nTip = 20)
  tc <- TransferConsensus(trees, greedy = "best")
  mr <- Consensus(trees, p = 0.5)

  # Transfer consensus should be at least as resolved
  expect_gte(NSplits(tc), NSplits(mr))
})

test_that("Both greedy strategies produce valid trees", {
  trees <- as.phylo(1:15, nTip = 12)

  tc_first <- TransferConsensus(trees, greedy = "first")
  tc_best  <- TransferConsensus(trees, greedy = "best")

  expect_s3_class(tc_first, "phylo")
  expect_s3_class(tc_best, "phylo")
  expect_equal(sort(TipLabels(tc_first)), sort(TipLabels(trees[[1]])))
  expect_equal(sort(TipLabels(tc_best)), sort(TipLabels(trees[[1]])))
})

test_that("scale = FALSE (unscaled) works", {
  trees <- as.phylo(1:10, nTip = 10)

  tc_scaled <- TransferConsensus(trees, scale = TRUE)
  tc_unscaled <- TransferConsensus(trees, scale = FALSE)

  expect_s3_class(tc_scaled, "phylo")
  expect_s3_class(tc_unscaled, "phylo")
})

test_that("init = 'majority' works", {
  trees <- as.phylo(0:19, nTip = 15)

  tc_empty <- TransferConsensus(trees, init = "empty")
  tc_maj   <- TransferConsensus(trees, init = "majority")

  expect_s3_class(tc_empty, "phylo")
  expect_s3_class(tc_maj, "phylo")
  # Both should produce reasonable trees
  expect_gte(NSplits(tc_empty), 1)
  expect_gte(NSplits(tc_maj), 1)
})

test_that("Error on bad input", {
  expect_error(TransferConsensus(list(rtree(5))), "multiPhylo")
  expect_error(TransferConsensus(structure(list(rtree(5)), class = "multiPhylo")),
               "at least 2")
})

test_that("Two-tree consensus returns a valid tree", {
  trees <- as.phylo(1:2, nTip = 8)
  tc <- TransferConsensus(trees)
  expect_s3_class(tc, "phylo")
  expect_equal(length(TipLabels(tc)), 8)
})
