test_that("Hierarchical Mutual Information", {
  skip_if_not_installed("TreeTools")
  library("TreeTools", quietly = TRUE)
  
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
  expect_error(HierarchicalMutualInfo(tree1), "tree2 must be provided")
  
  # Test with different tip numbers (should error)
  tree_small <- BalancedTree(6)
  expect_error(HierarchicalMutualInfo(tree1, tree_small))
  
  # Test reportMatching
  hmi_with_matching <- HierarchicalMutualInfo(tree1, tree2, reportMatching = TRUE)
  expect_true(is.numeric(hmi_with_matching))
  expect_true("matching" %in% names(attributes(hmi_with_matching)))
  
  # Test expected value for bal6 vs pec6 (should be approximately 0.24)
  bal6 <- BalancedTree(6)
  pec6 <- PectinateTree(6)
  hmi_bal_pec <- HierarchicalMutualInfo(bal6, pec6)
  
  # The expected value is 0.24 based on Python reference implementation
  expect_equal(hmi_bal_pec, 0.24, tolerance = 0.02)
}))

test_that("HMI helper functions", {
  skip_if_not_installed("TreeTools")
  library("TreeTools", quietly = TRUE)
  
  tree <- BalancedTree(8)
  
  # Test hierarchical partition building
  partition <- .PhyloToHierarchicalPartition(tree)
  
  expect_true(is.list(partition))
  
  # Test HMI recursive calculation
  tree2 <- PectinateTree(8)
  partition2 <- .PhyloToHierarchicalPartition(tree2)
  
  result <- .CalculateHMIRecursive(partition, partition2)
  expect_true(is.list(result))
  expect_true("n_ts" %in% names(result))
  expect_true("I_ts" %in% names(result))
  expect_true(is.numeric(result$n_ts))
  expect_true(is.numeric(result$I_ts))
  expect_true(result$n_ts >= 0)
  expect_true(result$I_ts >= 0)
}))

test_that("HMI comparison with standard mutual information", {
  library("TreeTools", quietly = TRUE)
  
  tree1 <- BalancedTree(8)
  tree2 <- PectinateTree(8)
  
  # Compare HMI with some basic principles
  hmi <- HierarchicalMutualInfo(tree1, tree2)
  
  # Both should be positive for different trees
  expect_true(hmi >= 0)
  
  # Test with identical trees
  hmi_identical <- HierarchicalMutualInfo(tree1, tree1)
  
  expect_true(hmi_identical >= 0)
  expect_true(is.numeric(hmi))
  expect_true(is.numeric(hmi_identical))
})

test_that("HMI with list inputs", {
  library("TreeTools", quietly = TRUE)
  
  trees <- list(
    BalancedTree(8),
    PectinateTree(8),
    RandomTree(8, 1)
  )
  
  # Test with list input
  hmi_result <- HierarchicalMutualInfo(trees)
  
  expect_true(inherits(hmi_result, "dist"))
  expect_equal(length(hmi_result), 3) # 3 pairwise distances for 3 trees
  
  # Convert to full matrix to test properties
  hmi_matrix <- as.matrix(hmi_result)
  expect_equal(dim(hmi_matrix), c(3, 3))
  
  # Matrix should be symmetric
  expect_equal(hmi_matrix[1, 2], hmi_matrix[2, 1], tolerance = 1e-10)
  expect_equal(hmi_matrix[1, 3], hmi_matrix[3, 1], tolerance = 1e-10)
  expect_equal(hmi_matrix[2, 3], hmi_matrix[3, 2], tolerance = 1e-10)
  
  # Diagonal should be zero (distance from tree to itself in distance matrix)
  expect_equal(hmi_matrix[1, 1], 0, tolerance = 1e-10)
  expect_equal(hmi_matrix[2, 2], 0, tolerance = 1e-10)
  expect_equal(hmi_matrix[3, 3], 0, tolerance = 1e-10)
})