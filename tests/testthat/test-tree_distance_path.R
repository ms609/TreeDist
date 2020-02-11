context('tree_distance_path.R')

test_that("path.dist doesn't crash", {
  library("TreeTools")
  expect_equal(c(5.66, 6, 6, 6.32, 6.32, 5.74),
               PathDist(as.phylo(0:5, 6), BalancedTree(6)),
               tolerance = 2)
})
