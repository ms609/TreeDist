context('tree_information.R')

library('TreeTools')

test_that('ClusteringInfo method works', {
  trees <- list(BalancedTree(8), PectinateTree(8))
  expect_equal(vapply(trees, ClusteringInfo, 0),
               ClusteringInfo(structure(trees, class = 'multiPhylo')))
  expect_equal(vapply(trees, ClusteringInfo, 0),
               ClusteringInfo(trees))
})
  