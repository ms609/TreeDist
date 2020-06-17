context('utilities.R')

test_that('SplitsInBinaryTree() works', {
  Tree <- ape::as.phylo
  expect_identical(5L, SplitsInBinaryTree(8))
  expect_identical(5L, SplitsInBinaryTree(Tree(0, 8)))
  expect_identical(rep(5L, 4), SplitsInBinaryTree(Tree(0:3, 8)))
})
