context('tree_distance_spr.R')

test_that("SPR.dist called safely", {
  library("TreeTools")
  PhangornSPR <- phangorn::SPR.dist
  expect_equal(PhangornSPR(structure(lapply(as.phylo(0:5, 6), Postorder),
                                     class = 'multiPhylo'),
                           Postorder(BalancedTree(6))),
               SPRDist(as.phylo(0:5, 6), BalancedTree(6)))
  expect_equivalent(SPRDist(BalancedTree(6), PectinateTree(6)),
                    SPRDist(list(BalancedTree(6), PectinateTree(6)))[1])
  
})
