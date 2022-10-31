library("TreeTools", quietly = TRUE)
bal8 <- BalancedTree(8)
pec8 <- PectinateTree(8)
bal8BG <- DropTip(bal8, c(1, 8))
pec8BG <- DropTip(pec8, c(1, 8))
bal8CH <- DropTip(bal8, 1:2)
pec8CH <- DropTip(pec8, 1:2)

test_that("Non-identical tips handled okay", {
  fullDist <- TreeDistance(bal8, pec8)
  expect_equal(TreeDistance(bal8, bal8BG), 0)
  expect_equal(MutualClusteringInfo(bal8, bal8BG), ClusteringEntropy(bal8BG))
  expect_equal(TreeDistance(bal8BG, bal8), 0)
  expect_equal(MutualClusteringInfo(bal8BG, bal8), ClusteringEntropy(bal8BG))
  expect_equal(TreeDistance(bal8BG, bal8CH), 0)
  expect_equal(MutualClusteringInfo(bal8BG, bal8CH),
               ClusteringEntropy(DropTip(bal8, c(1, 2, 8))))
  expect_equal(TreeDistance(bal8, pec8BG), TreeDistance(bal8BG, pec8BG))
  expect_equal(TreeDistance(bal8BG, pec8), TreeDistance(bal8BG, pec8BG))
  expect_equal(MutualClusteringInfo(bal8, pec8BG),
               MutualClusteringInfo(bal8BG, pec8BG))
  expect_equal(TreeDistance(BalancedTree(1:5), BalancedTree(6:9)), NaN)
  expect_equal(MutualClusteringInfo(BalancedTree(1:5), BalancedTree(6:9)), 0)
})
