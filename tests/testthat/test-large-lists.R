# Test to verify that the fix for large lists works
# This test focuses on the core fix without depending on external packages

test_that("robinson_foulds_all_pairs stack reset fix", {
  # Skip if TreeTools is not available 
  skip_if_not_installed("TreeTools")
  
  # Test basic functionality with smaller lists first
  set.seed(0)
  trees_small <- lapply(rep(6, 10), TreeTools::RandomTree, root = TRUE)
  expect_no_error(RobinsonFoulds(trees_small))
  
  # Test the specific case that was failing (128 vs 129 trees)
  # Use smaller trees to speed up testing
  trees_128 <- lapply(rep(6, 128), TreeTools::RandomTree, root = TRUE)
  expect_no_error(RobinsonFoulds(trees_128))
  
  trees_129 <- lapply(rep(6, 129), TreeTools::RandomTree, root = TRUE)
  expect_no_error(RobinsonFoulds(trees_129))
  
  # Test even larger list to ensure robustness
  trees_150 <- lapply(rep(6, 150), TreeTools::RandomTree, root = TRUE)
  expect_no_error(RobinsonFoulds(trees_150))
})