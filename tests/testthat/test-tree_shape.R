context("tree_shape.R")

test_that("Errors are handled", {
  expect_error(RootedTreeWithShape(as.integer64(-1)))
  expect_error(UnrootedTreeWithShape(31, 31))
})