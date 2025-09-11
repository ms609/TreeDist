library("TreeTools", quietly = TRUE)

test_that("Hierarchical Mutual Information", {
  
  # Create test trees
  tree1 <- BalancedTree(8)
  tree2 <- PectinateTree(8)
  tree3 <- BalancedTree(8) # Identical to tree1
  
  # Test basic functionality
  expect_no_error(HierarchicalMutualInfo(tree1, tree2))
  
  # Test that HMI is numeric and non-negative
  hmi <- HierarchicalMutualInfo(tree1, tree2)
  expect_true(is.numeric(hmi))
  expect_true(hmi >= 0)
  
  # Test normalization
  hmi_norm <- HierarchicalMutualInfo(tree1, tree2, normalize = TRUE)
  expect_true(is.numeric(hmi_norm))
  expect_true(hmi_norm >= 0)
  expect_true(hmi_norm <= 1)
  
  # Test symmetry
  hmi12 <- HierarchicalMutualInfo(tree1, tree2)
  hmi21 <- HierarchicalMutualInfo(tree2, tree1)
  expect_equal(hmi12, hmi21, tolerance = 1e-10)
  
  # Test identity property - normalized self-comparison should be 1
  hmi_self1_norm <- HierarchicalMutualInfo(tree1, tree1, normalize = TRUE)
  hmi_self2_norm <- HierarchicalMutualInfo(tree2, tree2, normalize = TRUE)
  
  expect_equal(hmi_self1_norm, 1, tolerance = 1e-10)
  expect_equal(hmi_self2_norm, 1, tolerance = 1e-10)
  
  # Test normalized identity
  hmi_self_norm <- HierarchicalMutualInfo(tree1, tree1, normalize = TRUE)
  expect_equal(hmi_self_norm, 1, tolerance = 1e-10)
  
  # Test error handling
  expect_equal(HierarchicalMutualInfo(tree1), SelfHMI(tree1))
  expect_warning(expect_equal(HierarchicalMutualInfo(tree1, norm = TRUE), 1),
                 "tree2")
  
  # Test with different tip numbers (should error)
  tree_small <- BalancedTree(6)
  expect_error(HierarchicalMutualInfo(tree1, tree_small),
               "number of leaves")
})

test_that("HMI edge cases", {
  bal9 <- BalancedTree(9)
  bal9b <- BalancedTree(paste0("t", c(3:1, 7:9, 6:4)))
  
  expect_lt(HMI(bal9, bal9b), HMI(bal9, bal9))
  
  expect_lt(HMI(bal9, bal9b, normalize = TRUE), 0.05)
  
  expect_equal(AHMI(StarTree(6), BalancedTree(6))[[1]], 0)
  expect_equal(AHMI(StarTree(2), BalancedTree(2)), structure(0, sem = NaN))
})

