test_that("HierarchicalMutualInfoDist basic functionality", {
  library("TreeTools", quietly = TRUE)
  
  # Create simple test trees
  tree1 <- ape::read.tree(text = "((a,b),(c,d));")
  tree2 <- ape::read.tree(text = "((a,c),(b,d));")
  tree3 <- ape::read.tree(text = "((a,b),(c,d));") # Same as tree1
  
  # Test basic functionality
  expect_is(HierarchicalMutualInfoDist(tree1, tree2), "numeric")
  expect_length(HierarchicalMutualInfoDist(tree1, tree2), 1)
  
  # Distance should be non-negative
  expect_gte(HierarchicalMutualInfoDist(tree1, tree2), 0)
  
  # Distance of a tree with itself should be 0
  expect_equal(HierarchicalMutualInfoDist(tree1, tree1), 0, tolerance = 1e-6)
  expect_equal(HierarchicalMutualInfoDist(tree1, tree3), 0, tolerance = 1e-6)
  
  # Distance should be symmetric
  expect_equal(HierarchicalMutualInfoDist(tree1, tree2), 
               HierarchicalMutualInfoDist(tree2, tree1), tolerance = 1e-6)
})

test_that("HierarchicalMutualInfoDist with different tree sizes", {
  library("TreeTools", quietly = TRUE)
  
  # Test with balanced trees of different sizes
  tree4 <- BalancedTree(4)
  tree8 <- BalancedTree(8)
  
  # Should handle different sized trees by reducing to common tips
  # (This behavior should match other distance functions)
  expect_is(HierarchicalMutualInfoDist(tree4, tree8), "numeric")
  
  # Test with star trees
  star4 <- StarTree(4)
  star8 <- StarTree(8)
  
  expect_is(HierarchicalMutualInfoDist(star4, star8), "numeric")
})

test_that("HierarchicalMutualInfoDist normalization", {
  library("TreeTools", quietly = TRUE)
  
  tree1 <- BalancedTree(6)
  tree2 <- PectinateTree(6)
  
  # Test normalization
  unnormalized <- HierarchicalMutualInfoDist(tree1, tree2, normalize = FALSE)
  normalized <- HierarchicalMutualInfoDist(tree1, tree2, normalize = TRUE)
  
  expect_is(unnormalized, "numeric")
  expect_is(normalized, "numeric")
  
  # Normalized should be between 0 and 1
  expect_gte(normalized, 0)
  expect_lte(normalized, 1)
  
  # Normalized should be different from unnormalized (unless distance is 0)
  if (unnormalized > 0) {
    expect_ne(unnormalized, normalized)
  }
})

test_that("HierarchicalMutualInfoDist with lists of trees", {
  library("TreeTools", quietly = TRUE)
  
  # Create list of trees
  trees <- list(
    BalancedTree(6),
    PectinateTree(6),
    StarTree(6)
  )
  
  # Test pairwise distances
  distances <- HierarchicalMutualInfoDist(trees)
  
  expect_is(distances, "matrix")
  expect_equal(nrow(distances), 3)
  expect_equal(ncol(distances), 3)
  
  # Diagonal should be 0 (distance of tree with itself)
  expect_equal(diag(distances), rep(0, 3), tolerance = 1e-6)
  
  # Matrix should be symmetric
  expect_equal(distances[1, 2], distances[2, 1], tolerance = 1e-6)
  expect_equal(distances[1, 3], distances[3, 1], tolerance = 1e-6)
  expect_equal(distances[2, 3], distances[3, 2], tolerance = 1e-6)
})

test_that("HierarchicalMutualInfoDist edge cases", {
  library("TreeTools", quietly = TRUE)
  
  # Test with minimum tree (3 tips)
  tree3 <- ape::read.tree(text = "(a,(b,c));")
  expect_is(HierarchicalMutualInfoDist(tree3, tree3), "numeric")
  expect_equal(HierarchicalMutualInfoDist(tree3, tree3), 0, tolerance = 1e-6)
  
  # Test with star tree (no internal structure)
  star <- StarTree(5)
  balanced <- BalancedTree(5)
  
  expect_is(HierarchicalMutualInfoDist(star, balanced), "numeric")
  expect_gte(HierarchicalMutualInfoDist(star, balanced), 0)
})

test_that("HierarchicalMutualInfoDist matches expected values for simple cases", {
  library("TreeTools", quietly = TRUE)
  
  # Test case 1: Two completely different 4-tip trees
  tree1 <- ape::read.tree(text = "((a,b),(c,d));")
  tree2 <- ape::read.tree(text = "((a,c),(b,d));")
  
  dist1 <- HierarchicalMutualInfoDist(tree1, tree2)
  expect_is(dist1, "numeric")
  expect_gt(dist1, 0) # Should be positive for different trees
  
  # Test case 2: Star tree vs balanced tree
  star4 <- StarTree(4)
  balanced4 <- BalancedTree(4)
  
  dist2 <- HierarchicalMutualInfoDist(star4, balanced4)
  expect_is(dist2, "numeric")
  expect_gt(dist2, 0) # Should be positive
  
  # The star tree should be maximally different from balanced tree
  # This tests that our algorithm captures hierarchical differences
  expect_gt(dist2, dist1) # Star vs balanced should be more different
})

test_that("HierarchicalMutualInfoDist consistency with other distance functions", {
  library("TreeTools", quietly = TRUE)
  
  # Compare behavior with ClusteringInfoDistance on same trees
  tree1 <- BalancedTree(6)
  tree2 <- PectinateTree(6)
  
  hmi_dist <- HierarchicalMutualInfoDist(tree1, tree2)
  clustering_dist <- ClusteringInfoDistance(tree1, tree2)
  
  # Both should be non-negative
  expect_gte(hmi_dist, 0)
  expect_gte(clustering_dist, 0)
  
  # Both should give 0 for identical trees
  expect_equal(HierarchicalMutualInfoDist(tree1, tree1), 0, tolerance = 1e-6)
  expect_equal(ClusteringInfoDistance(tree1, tree1), 0, tolerance = 1e-6)
  
  # Both should be symmetric
  expect_equal(HierarchicalMutualInfoDist(tree1, tree2),
               HierarchicalMutualInfoDist(tree2, tree1), tolerance = 1e-6)
  expect_equal(ClusteringInfoDistance(tree1, tree2),
               ClusteringInfoDistance(tree2, tree1), tolerance = 1e-6)
})