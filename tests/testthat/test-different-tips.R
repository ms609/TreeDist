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
  
  fullList <- list(
    bal8 = bal8,
    pec8 = pec8,
    
    bal8BG = bal8BG,
    pec8BG = pec8BG,
    bal8CH = bal8CH,
    pec8CH = pec8CH,
    
    balAE = BalancedTree(5),
    balFI = BalancedTree(paste0("t", 6:9))
  )
  
  expect_equal(MutualClusteringInfo(bal8, fullList),
               unlist(vapply(fullList, MutualClusteringInfo, 1, bal8)))
  expect_equal(MutualClusteringInfo(bal8, fullList, normalize = TRUE),
               unlist(vapply(fullList, MutualClusteringInfo, 1, bal8,
                             normalize = TRUE)))
  expect_equal(TreeDistance(bal8, fullList),
               unlist(lapply(fullList, TreeDistance, bal8)))
  expect_equal(TreeDistance(fullList, bal8), TreeDistance(bal8, fullList))
  expect_equal(MutualClusteringInfo(fullList, bal8),
               MutualClusteringInfo(bal8, fullList))
  expect_equal(MutualClusteringInfo(fullList, fullList),
               vapply(fullList, function(t1)
                 vapply(fullList, MutualClusteringInfo, double(1), t1),
                 double(length(fullList)))
               )
  expect_equal(MutualClusteringInfo(fullList),
               vapply(fullList, MutualClusteringInfo,
                      double(length(fullList)), fullList))
  expect_equal(MutualClusteringInfo(fullList, normalize = TRUE),
               vapply(fullList, MutualClusteringInfo,
                      double(length(fullList)), fullList, normalize = TRUE))
  expect_equal(TreeDistance(fullList),
               vapply(fullList, TreeDistance, double(length(fullList)),
                      fullList)))
})
