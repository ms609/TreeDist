library(ape)

context("tree rearrangements")
test_that("replacement reorder functions work correctly", {
  ## Tree
  tree <- read.tree(text = "((((((1,2),3),4),5),6),(7,(8,(9,(10,(11,12))))));")
  expect_equal(Cladewise(tree), reorder(tree, 'cladewise'))
  expect_equal(Pruningwise(tree), reorder(tree, 'pruningwise'))
  expect_equal(Postorder(tree), reorder(tree, 'postorder'))
})
