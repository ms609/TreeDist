test_that("Hierarchical Mutual Information", {
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
  
  # Test with splits objects
  splits1 <- as.Splits(tree1)
  splits2 <- as.Splits(tree2)
  
  hmi_splits <- HierarchicalMutualInfoSplits(splits1, splits2)
  expect_equal(hmi, hmi_splits, tolerance = 1e-10)
  
  # Test error handling
  expect_error(HierarchicalMutualInfo(tree1), "tree2 must be provided")
  
  # Test with different tip numbers (should error)
  tree_small <- BalancedTree(6)
  expect_error(HierarchicalMutualInfo(tree1, tree_small))
  
  # Test reportMatching
  hmi_with_matching <- HierarchicalMutualInfo(tree1, tree2, reportMatching = TRUE)
  expect_true(is.numeric(hmi_with_matching))
  expect_true("matching" %in% names(attributes(hmi_with_matching)))
})

test_that("HMI helper functions", {
  library("TreeTools", quietly = TRUE)
  
  tree <- BalancedTree(8)
  splits <- as.Splits(tree)
  
  # Test hierarchical weights calculation
  weights <- .CalculateHierarchicalWeights(splits, 8)
  
  expect_true(is.numeric(weights))
  expect_equal(length(weights), length(splits))
  expect_true(all(weights >= 0))
  expect_equal(sum(weights), 1, tolerance = 1e-10) # Should be normalized
  
  # Test with empty splits - just create an empty matrix  
  empty_splits <- matrix(raw(0), 0, 0)
  empty_weights <- .CalculateHierarchicalWeights(empty_splits, 0)
  expect_equal(length(empty_weights), 0)
  
  # Test weighted mutual information calculation
  tree2 <- PectinateTree(8)
  splits2 <- as.Splits(tree2)
  weights2 <- .CalculateHierarchicalWeights(splits2, 8)
  
  wmi <- .CalculateWeightedMutualInfo(splits, splits2, weights, weights2, 8)
  expect_true(is.numeric(wmi))
  expect_true(wmi >= 0)
  
  # Test maximum HMI calculation - max should be based on self-comparison
  max_hmi1 <- .MaxHierarchicalMutualInfo(splits, splits, weights, weights)
  max_hmi2 <- .MaxHierarchicalMutualInfo(splits2, splits2, weights2, weights2)
  expect_true(is.numeric(max_hmi1))
  expect_true(max_hmi1 >= 0)
  expect_true(is.numeric(max_hmi2))
  expect_true(max_hmi2 >= 0)
})

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