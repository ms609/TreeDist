library('TreeTools')

test_that("SplitwiseInfo() handles probabilities", {
  Tree <- function (txt) ape::read.tree(text = txt)
  tree <- Tree('((a, b)60, (c, d));')
  expect_error(SplitwiseInfo(tree, TRUE))
  expect_gt(SplitwiseInfo(tree), SplitwiseInfo(tree, 100))
  expect_equal(0, SplitwiseInfo(tree, 60 * 3),
               tolerance = sqrt(.Machine$double.eps))
  
  p <- 1.0 * c(1, 0, 0) + 0.0 * c(0, 1/2, 1/2)
  expect_equal(log2(3), SplitwiseInfo(tree))
  
  p <- 0.6 * c(1, 0, 0) + 0.4 * c(0, 1/2, 1/2)
  expect_equal(sum(p), 1)
  expectation <- log2(3) + sum(p * log2(p))
  expect_equal(expectation, SplitwiseInfo(tree, 100))
  
  
  treeP <- Tree('((a, b)0.60, (c, d));')
  expect_equal(SplitwiseInfo(tree, 100), SplitwiseInfo(treeP, TRUE))
  
  expect_equal(SplitwiseInfo(Tree('(a, b, (c, (d, e)0.8)0.75);'), TRUE),
               SplitwiseInfo(Tree('(a, b, (c, d, e)0.75);'), TRUE) +
                 SplitwiseInfo(Tree('(a, b, c, (d, e)0.8);'), TRUE))
  expect_equal(SplitwiseInfo(Tree('(a, b, (c, (d, e)0.8)0.75);'), TRUE),
               SplitwiseInfo(Tree('(a, b, (c, (d, e)));'), c(0.75, 0.8)))
  expect_equal(SplitwiseInfo(Tree('(a, b, (c, (d, e)0.8));'), TRUE),
               SplitwiseInfo(Tree('(a, b, (c, (d, e)));'), c(1, 0.8)))
  expect_equal(SplitwiseInfo(Tree('(a, b, (c, (d, e)));')),
               SplitwiseInfo(Tree('(a, b, (c, (d, e)));'), TRUE))
  expect_equal(SplitwiseInfo(Tree('(a, b, (c, (d, e)));'), TRUE),
               SplitwiseInfo(Tree('(a, b, (c, (d, e)));'), c(1, 1)))  
  
})

test_that('ClusteringInfo() method works', {
  trees <- list(BalancedTree(8), PectinateTree(8))
  expect_equal(vapply(trees, ClusteringInfo, 0),
               ClusteringInfo(structure(trees, class = 'multiPhylo')))
  expect_equal(vapply(trees, ClusteringInfo, 0),
               ClusteringInfo(trees))
})
  
test_that("ConsensusInfo() is robust", {
  trees <- list(ape::read.tree(text = '(a, (b, (c, (d, (e, X)))));'),
                ape::read.tree(text = '((a, X), (b, (c, (d, e))));'))
  expect_equal(0, ConsensusInfo(trees, 'cl'))
})
