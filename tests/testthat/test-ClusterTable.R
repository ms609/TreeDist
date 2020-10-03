context("ClusterTable.R")
library("TreeTools", quietly = TRUE, warn.conflicts = FALSE)

test_that("ClusterTable class behaves", {
  tree <- RootTree(BalancedTree(6), 1)
  ct <- as.ClusterTable(tree)
  expect_equal(matrix(c(0, 2, rep(1, 4), 0, 3, 3:6), 6),
               as.matrix.ClusterTable(ct))
})

test_that("Attributes are correct", {
  t6 <- as.ClusterTable(BalancedTree(6))
  t7 <- as.ClusterTable(PectinateTree(7))
  t8 <- as.ClusterTable(BalancedTree(8))
  s8 <- StarTree(8)
  expect_equal(3, NSplits(t6))
  expect_equal(4:5, NSplits(list(t7, t8)))
  
  expect_equal(6, NTip(t6))
  expect_equal(7:8, NTip(list(t7, t8)))
  
  #TODO test TipLabels, SplitsInBalancedTree
})