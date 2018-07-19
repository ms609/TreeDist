context("Ape functions")
test_that("File time is read correctly", {
  expect_equal('2018-07-18 13:47:46', ApeTime('test-ape-tree.nex', 'string'))
})