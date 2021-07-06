library('TreeTools')

test_that("SplitwiseInfo() handles probabilities", {
  Tree <- function (txt) ape::read.tree(text = txt)
  tree <- Tree('((a, b)60, (c, d));')
  Test <- function (Expect, tree, p, ...) {
    Expect(..., SplitwiseInfo(tree, p))
    Expect(..., ClusteringInfo(tree, p))
    Expect(..., ClusteringEntropy(tree, p))
  }
  Clust <- function (tree, ...) {
    expect_equal(ClusteringInfo(tree, ..., sum = TRUE),
                 sum(ClusteringInfo(tree, ..., sum = FALSE)))
    expect_equal(ClusteringEntropy(tree, ..., sum = TRUE),
                 sum(ClusteringEntropy(tree, ..., sum = FALSE)))
    expect_equal(ClusteringInfo(tree, ...) / NTip(tree),
                 ClusteringEntropy(tree, ...))
  }
  Test(expect_error, tree, TRUE)
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
  expect_equal(Clust(Tree('(a, b, (c, (d, e)));'), TRUE),
               Clust(Tree('(a, b, (c, (d, e)));'), c(1, 1)))
  
  expect_equal(SplitwiseInfo(Tree('(a, b, (c, (d, e)0.8));'), TRUE),
               sum(SplitwiseInfo(Tree('(a, b, (c, (d, e)0.8));'), TRUE, FALSE)))
})

test_that("SplitwiseInfo() can't be improved by dropping resolved tip", {
  # Test arguably redundant, but a useful reminder of a closed optimization
  # possibility
  b8With <- as.Splits(c(T, T, T, T, F, F, F, F))
  b8Without <- as.Splits(c(T, T, T, F, F, F, F))
  i8With <- as.Splits(c(T, T, T, T, T, T, F, F))
  i8Without <- as.Splits(c(T, T, T, T, T, F, F))
  expect_lt(SplitwiseInfo(b8With, p = 0.5), SplitwiseInfo(b8Without))
  expect_lt(SplitwiseInfo(i8With, p = 0.5), SplitwiseInfo(i8Without))

})

test_that('ClusteringInfo() method works', {
  trees <- list(BalancedTree(8), PectinateTree(8))
  expect_equal(vapply(trees, ClusteringInfo, 0),
               ClusteringInfo(structure(trees, class = 'multiPhylo')))
  expect_equal(vapply(trees, ClusteringInfo, 0),
               ClusteringInfo(trees))
})
  