context("mast.cpp")

test_that('MAST works', {
  tree1 <- BalancedTree(8)
  tree2 <- PectinateTree(8)
  MAST(tree1, tree2)
})