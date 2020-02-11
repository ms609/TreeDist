context('tree_distance_path.R')

test_that("path.dist called safely", {
  library("TreeTools")
  expect_equal(c(5.66, 6, 6, 6.32, 6.32, 5.74),
               PathDist(as.phylo(0:5, 6), BalancedTree(6)),
               tolerance = 2)
  
  expect_equivalent(PathDist(BalancedTree(6), PectinateTree(6)),
                    PathDist(list(BalancedTree(6), PectinateTree(6)))[1])
})
