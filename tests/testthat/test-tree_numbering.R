context("tree_numbering.R")
test_that("replacement reorder functions work correctly", {
  ## Tree
  tree <- ape::read.tree(text = 
                           "((((((1,2),3),4),5),6),(7,(8,(9,(10,(11,12))))));")
  expect_equal(Cladewise(tree), ape::reorder.phylo(tree, 'cladewise'))
  expect_equal(Pruningwise(tree), ape::reorder.phylo(tree, 'pruningwise'))
  expect_equal(Postorder(tree), ape::reorder.phylo(tree, 'postorder'))
})
