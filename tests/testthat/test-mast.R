context("mast.cpp")
library('TreeTools')

test_that('MAST works', {
  tree1 <- BalancedTree(8L)
  tree2 <- PectinateTree(8L)
  expect_equal(8L, MASTSize(tree1, tree1, rooted = TRUE))
  expect_equal(8L, MASTSize(tree1, tree1, rooted = FALSE))
  expect_equal(4L, MASTSize(tree1, tree2, rooted = TRUE))
  expect_equal(6L, MASTSize(tree1, tree2, rooted = FALSE))
  
  expect_equal(MASTSize(BalancedTree(7), as.phylo(0:3, 7)),
               MASTSize(as.phylo(0:3, 7), BalancedTree(7)))
  
  expect_equal(MASTSize(list(BalancedTree(7), PectinateTree(7)), as.phylo(0:3, 7))[1, ],
               MASTSize(as.phylo(0:3, 7), BalancedTree(7)))
  
  expect_equal(MASTInfo(BalancedTree(7), as.phylo(0:3, 7)),
               MASTInfo(as.phylo(0:3, 7), BalancedTree(7)))
  
  expect_equal(MASTInfo(list(BalancedTree(7), PectinateTree(7)), as.phylo(0:3, 7))[1, ],
               MASTInfo(as.phylo(0:3, 7), BalancedTree(7)))
})
