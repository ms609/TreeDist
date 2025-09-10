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

test_that("HMI with multiple trees", {
  # Test with list of trees
  trees <- list(bal = BalancedTree(6), pec = PectinateTree(6), star = StarTree(6))
  
  # Test pairwise distances  
  dist_result <- HierarchicalMutualInfo(trees)
  expect_s3_class(dist_result, "dist")
  expect_equal(attr(dist_result, "Size"), 3)
  expect_equal(attr(dist_result, "Labels"), names(trees))
  
  # Test single tree vs list
  single_vs_list <- HierarchicalMutualInfo(trees[[1]], trees[-1])
  expect_length(single_vs_list, 2)
  expect_true(all(single_vs_list >= 0))
  
  # Test list vs list  
  trees2 <- list(BalancedTree(6), PectinateTree(6))
  list_vs_list <- HierarchicalMutualInfo(trees[1:2], trees2)
  expect_equal(dim(list_vs_list), c(2, 2))
  
  # Test normalization with multiple trees
  dist_norm <- HierarchicalMutualInfo(trees, normalize = TRUE)
  expect_s3_class(dist_norm, "dist")
  expect_true(all(dist_norm >= 0))
  expect_true(all(dist_norm <= 1))
})

test_that("HMI with larger trees", {
  # Test with larger trees as requested by @ms609
  bal8 <- BalancedTree(8) 
  pec8 <- PectinateTree(8)
  star8 <- StarTree(8)
  
  expect_no_error(HierarchicalMutualInfo(bal8, pec8))
  expect_no_error(HierarchicalMutualInfo(bal8, star8))
  expect_no_error(HierarchicalMutualInfo(pec8, star8))
  
  # Test with random trees
  set.seed(123)
  rnd8 <- RandomTree(8, root = TRUE)
  
  expect_no_error(HierarchicalMutualInfo(bal8, rnd8))
  expect_no_error(HierarchicalMutualInfo(rnd8, pec8))
  
  # Values should be reasonable
  hmi_bal_pec <- HierarchicalMutualInfo(bal8, pec8)
  hmi_bal_rnd <- HierarchicalMutualInfo(bal8, rnd8)
  hmi_rnd_pec <- HierarchicalMutualInfo(rnd8, pec8)
  
  expect_true(hmi_bal_pec > 0)
  expect_true(hmi_bal_rnd > 0)
  expect_true(hmi_rnd_pec > 0)
})

test_that("HMI edge cases", {
  bal9 <- BalancedTree(9)
  bal9b <- BalancedTree(paste0("t", c(3:1, 7:9, 6:4)))
  
  expect_lt(HMI(bal9, bal9b), HMI(bal9, bal9))
  
  expect_lt(HMI(bal9, bal9b, normalize = TRUE), 0.05)
  
  expect_equal(AHMI(StarTree(6), BalancedTree(6))[[1]], 0)
  expect_equal(AHMI(StarTree(2), BalancedTree(2)), structure(NaN, sem = NaN))
})

test_that("EHMI and AHMI functions", {
  tree1 <- BalancedTree(6)
  tree2 <- PectinateTree(6)
  
  # Test EHMI
  ehmi_result <- EHMI(tree1, tree2, tolerance = 0.1, minResample = 10)
  expect_true(is.numeric(ehmi_result))
  expect_true(length(ehmi_result) == 1)
  expect_true(!is.null(attr(ehmi_result, "var")))
  expect_true(!is.null(attr(ehmi_result, "sem")))
  
  # Test AHMI
  ahmi_result <- AHMI(tree1, tree2, tolerance = 0.1, minResample = 10)
  expect_true(is.numeric(ahmi_result))
  expect_true(length(ahmi_result) == 1)
  expect_true(!is.null(attr(ahmi_result, "sem")))
})

test_that("SelfHMI function", {
  tree <- BalancedTree(6)
  
  # Test SelfHMI
  self_hmi <- SelfHMI(tree)
  expect_true(is.numeric(self_hmi))
  expect_true(self_hmi > 0)
  
  # Should be same as HMI with itself
  expect_equal(self_hmi, HierarchicalMutualInfo(tree), tolerance = 1e-10)
})

test_that("Error handling", {
  tree1 <- BalancedTree(6)
  tree2 <- BalancedTree(8)
  
  # Different number of tips should error
  expect_error(HierarchicalMutualInfo(tree1, tree2))
  
  # Invalid normalize parameter
  expect_error(HierarchicalMutualInfo(tree1, tree1, normalize = "invalid"))
})

