context("mast.cpp")
library('TreeTools')

test_that('MAST works', {
  tree1 <- BalancedTree(8L)
  tree2 <- PectinateTree(8L)
  expect_equal(8L, MAST(tree1, tree1, rooted = TRUE))
  expect_equal(8L, MAST(tree1, tree1, rooted = FALSE))
  expect_equal(4L, MAST(tree1, tree2, rooted = TRUE))
  expect_equal(6L, MAST(tree1, tree2, rooted = FALSE))
})