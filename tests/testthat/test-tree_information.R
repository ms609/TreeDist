library('TreeTools')

test_that("SplitwiseInfo() / ClusteringInfo() handle probabilities", {
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
  expect_gt(ClusteringInfo(tree), ClusteringInfo(tree, 100))
  expect_gt(ClusteringEntropy(tree), ClusteringEntropy(tree, 100))
  
  expect_equal(0, SplitwiseInfo(tree, 60 * 3),
               tolerance = sqrt(.Machine$double.eps))
  expect_equal(1 / 3, Clust(tree, 60 * 3))
  
  p <- 1.0 * c(1, 0, 0) + 0.0 * c(0, 1/2, 1/2)
  expect_equal(log2(3), SplitwiseInfo(tree))
  expect_equal(1, Clust(tree))
  
  p <- 0.6 * c(1, 0, 0) + 0.4 * c(0, 1/2, 1/2)
  expect_equal(sum(p), 1)
  expectation <- log2(3) + sum(p * log2(p))
  expect_equal(expectation, SplitwiseInfo(tree, 100))
  expect_equal(0.6, Clust(tree, 100))
  
  
  treeP <- Tree('((a, b)0.60, (c, d));')
  expect_equal(SplitwiseInfo(tree, 100), SplitwiseInfo(treeP, TRUE))
  expect_equal(Clust(tree, 100), Clust(treeP, TRUE))
  
  expect_equal(SplitwiseInfo(Tree('(a, b, (c, (d, e)0.8)0.75);'), TRUE),
               SplitwiseInfo(Tree('(a, b, (c, d, e)0.75);'), TRUE) +
                 SplitwiseInfo(Tree('(a, b, c, (d, e)0.8);'), TRUE))
  expect_equal(Clust(Tree('(a, b, (c, (d, e)0.8)0.75);'), TRUE),
               Clust(Tree('(a, b, (c, d, e)0.75);'), TRUE) +
                 Clust(Tree('(a, b, c, (d, e)0.8);'), TRUE))
  expect_equal(SplitwiseInfo(Tree('(a, b, (c, (d, e)0.8)0.75);'), TRUE),
               SplitwiseInfo(Tree('(a, b, (c, (d, e)));'), c(0.75, 0.8)))
  expect_equal(Clust(Tree('(a, b, (c, (d, e)0.8)0.75);'), TRUE),
               Clust(Tree('(a, b, (c, (d, e)));'), c(0.75, 0.8)))
  expect_equal(SplitwiseInfo(Tree('(a, b, (c, (d, e)0.8));'), TRUE),
               SplitwiseInfo(Tree('(a, b, (c, (d, e)));'), c(1, 0.8)))
  expect_equal(Clust(Tree('(a, b, (c, (d, e)0.8));'), TRUE),
               Clust(Tree('(a, b, (c, (d, e)));'), c(1, 0.8)))
  expect_equal(SplitwiseInfo(Tree('(a, b, (c, (d, e)));')),
               SplitwiseInfo(Tree('(a, b, (c, (d, e)));'), TRUE))
  expect_equal(Clust(Tree('(a, b, (c, (d, e)));')),
               Clust(Tree('(a, b, (c, (d, e)));'), TRUE))
  expect_equal(SplitwiseInfo(Tree('(a, b, (c, (d, e)));'), TRUE),
               SplitwiseInfo(Tree('(a, b, (c, (d, e)));'), c(1, 1)))  
  expect_equal(Clust(Tree('(a, b, (c, (d, e)));'), TRUE),
               Clust(Tree('(a, b, (c, (d, e)));'), c(1, 1)))  
})


test_that("SplitwiseInfo() / ClusteringInfo(sum = FALSE)", {
  splits <- as.Splits(BalancedTree(8))
  x <- SplitwiseInfo(BalancedTree(8), sum = FALSE)
  expect_equal(length(x), length(splits))
  expect_equal(names(x), names(splits))
  x <- ClusteringInfo(BalancedTree(8), sum = FALSE)
  expect_equal(length(x), length(splits))
  expect_equal(names(x), names(splits))
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
  
  balancedWithout <- BalancedTree(32)
  balancedWith <- AddTip(balancedWithout, 32)
  p <- double(30) + 1
  p[c(58, 62, 64, 65) - 35] <- 0.5
  expect_lt(SplitwiseInfo(balancedWith, p = p), SplitwiseInfo(balancedWithout))
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
