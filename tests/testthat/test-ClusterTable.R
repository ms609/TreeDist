context("ClusterTable.R")
library("TreeTools", quietly = TRUE)

test_that("ClusterTable class behaves", {
  tree <- BalancedTree(6)
  as.ClusterTable(tree)
  
})