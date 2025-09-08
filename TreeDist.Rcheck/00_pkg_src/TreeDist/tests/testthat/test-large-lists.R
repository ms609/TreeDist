test_that("RobinsonFoulds() supports large lists", {
  # Test that was failing before the fix
  set.seed(0)
  trees <- lapply(rep(8, 128), TreeTools::RandomTree, root = TRUE)
  expect_no_error(RobinsonFoulds(trees)) # Should work
  
  trees <- lapply(rep(8, 129), TreeTools::RandomTree, root = TRUE)
  expect_no_error(RobinsonFoulds(trees)) # Should work now
  
  # Test even larger list to ensure robustness
  trees <- lapply(rep(8, 150), TreeTools::RandomTree, root = TRUE)
  expect_no_error(RobinsonFoulds(trees))
})