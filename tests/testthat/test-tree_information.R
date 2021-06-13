library('TreeTools')

test_that("SplitwiseInfo() handles probabilities", {
  Tree <- function (txt) ape::read.tree(text = txt)
  tree <- Tree('((a, b)60, (c, d));')
  expect_error(SplitwiseInfo(tree, TRUE))
  expect_gt(SplitwiseInfo(tree), SplitwiseInfo(tree, 100))
  
  p <- 0.6 * c(1, 0, 0) + 0.4 * c(0, 1/2, 1/2)
  expect_equal(sum(p), 1)
  expectation <- -sum(p * log2(p))
  expect_equal(-expectation, SplitwiseInfo(tree, 100))
  
  
  treeP <- Tree('((a, b)0.60, (c, d));')
  expect_equal(SplitwiseInfo(tree, 100), SplitwiseInfo(treeP, TRUE))
  
  expect_equal(SplitwiseInfo(Tree('(a, b, (c, (d, e)0.8)0.75);'), TRUE),
               SplitwiseInfo(Tree('(a, b, (c, d, e)0.75);'), TRUE) +
                 SplitwiseInfo(Tree('(a, b, c, (d, e)0.8);'), TRUE))
               
               
  
})

test_that('ClusteringInfo() method works', {
  trees <- list(BalancedTree(8), PectinateTree(8))
  expect_equal(vapply(trees, ClusteringInfo, 0),
               ClusteringInfo(structure(trees, class = 'multiPhylo')))
  expect_equal(vapply(trees, ClusteringInfo, 0),
               ClusteringInfo(trees))
})
  