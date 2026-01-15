test_that("Islands() works", {
  library("TreeTools", quietly = TRUE)
  trees <- as.phylo(as.TreeNumber(BalancedTree(16)) + c(70:78, -(39:30), 80:99), 16)
  distances <- ClusteringInfoDist(trees)
  expect_equal(Islands(distances, 2.5), rep(c(1, 2, 1), c(9, 10, 20)))
  expect_equal(Islands(distances, 2.5, FALSE), rep(c(1, 10, 1), c(9, 10, 20)))
  
  
  trees <- as.phylo(c(0:10, 1000:1005, 2000:2019), 16)
  distances <- ClusteringInfoDist(trees)
  expect_equal(Islands(distances, 2.5, TRUE, smallest = 6),
               rep(c(1, 2, 3), c(11, 6, 20)))
  expect_equal(Islands(distances, 2.5, TRUE, smallest = 7),
               rep(c(1, NA, 2), c(11, 6, 20)))
  expect_equal(Islands(distances, 2.5, FALSE, smallest = 6),
               rep(c(1, 12, 18), c(11, 6, 20)))
  expect_equal(Islands(distances, 2.5, FALSE, smallest = 7),
               rep(c(1, NA, 18), c(11, 6, 20)))
})
