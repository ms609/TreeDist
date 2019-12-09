context("mast.cpp")
library('TreeTools')

test_that('MAST works', {
  tree1 <- BalancedTree(8)
  tree2 <- PectinateTree(8)
  MAST(tree1, tree2)
  expect_equal(nTip, cpp_mast(tree1$edge - 1L, tree1$edge - 1L, nTip))
  expect_equal(4, cpp_mast(tree1$edge - 1L, tree2$edge - 1L, nTip))
  
  expect_equal(4, MAST(tree1, tree2, rooted = TRUE))
  expect_equal(5, MAST(tree1, tree2, rooted = FALSE))
})