context('tree_distance_spr.R')

test_that("SPR.dist called safely", {
  library("TreeTools")
  expect_equal(c(1,1,1,1,2,2), SPRDist(as.phylo(0:5, 6), BalancedTree(6)))
  expect_equivalent(SPRDist(BalancedTree(6), PectinateTree(6)),
                    SPRDist(list(BalancedTree(6), PectinateTree(6)))[1])
  
})
