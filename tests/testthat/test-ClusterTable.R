context("ClusterTable.R")
library("TreeTools", quietly = TRUE, warn.conflicts = FALSE)

test_that("ClusterTable class behaves", {
  tree <- RootTree(BalancedTree(6), 1)
  ct <- as.ClusterTable(tree)
  expect_equal(matrix(c(0, 2:5, 1, rep(0, 5), rep(6, 5), rep(0, 4)), 10),
               as.matrix.ClusterTable(ct))
})