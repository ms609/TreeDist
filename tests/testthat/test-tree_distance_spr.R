context('tree_distance_spr.R')

test_that("SPR.dist doesn't crash", {
  library("TreeTools")
  expect_equal(c(1,1,1,1,2,2), SPRDist(as.phylo(0:5, 6), BalancedTree(6)))
})
