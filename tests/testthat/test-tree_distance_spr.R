context('tree_distance_spr.R')

test_that("SPR.dist called safely", {
  library("TreeTools")
  PhangornSPR <- phangorn::SPR.dist
  expect_equal(PhangornSPR(structure(lapply(as.phylo(0:5, 6), Postorder),
                                     class = 'multiPhylo'),
                           Postorder(BalancedTree(6))),
               SPRDist(as.phylo(0:5, 6), BalancedTree(6)))
  expect_equal(SPRDist(BalancedTree(6), as.phylo(0:5, 6)),
               SPRDist(as.phylo(0:5, 6), BalancedTree(6)))
  expect_equivalent(SPRDist(BalancedTree(6), PectinateTree(6)),
                    SPRDist(list(BalancedTree(6), PectinateTree(6)))[1])
  
  # https://github.com/KlausVigo/phangorn/issues/97
  tr1 <- structure(list(edge = structure(c(11L, 11L, 10L, 10L, 9L, 9L, 8L, 8L, 7L,
                                           7L, 2L, 6L, 5L, 11L, 4L, 10L, 3L, 9L, 
                                           1L, 8L), .Dim = c(10L, 2L)),
                        tip.label = c("t1", "t2", "t3", "t4", "t5", "t6"), 
                        Nnode = 5),
                   class = "phylo")
  
  tr2 <- structure(list(edge = structure(c(10L, 10L, 11L, 11L, 9L, 9L, 8L, 8L, 7L,
                                           7L, 2L, 6L, 3L, 4L, 5L, 10L, 11L, 9L, 
                                           1L, 8L), .Dim = c(10L, 2L)),
                        tip.label = c("t1", "t2", "t3", "t4", "t5", "t6"),
                        Nnode = 5),
                   class = "phylo")
  
  tr3 <- structure(list(edge = structure(c(9L, 9L, 11L, 11L, 10L, 10L, 8L, 8L, 7L, 
                                           7L, 1L, 2L, 4L, 5L, 6L, 11L, 3L, 9L, 
                                           8L, 10L), .Dim = c(10L, 2L)),
                        Nnode = 5L, 
                        tip.label = c("t1", "t2", "t3", "t4", "t5", "t6")),
                   class = "phylo")
  trs12 <- structure(list(tr1, tr2), class = 'multiPhylo')
  trs123 <- structure(list(tr1, tr2, tr3), class = 'multiPhylo')
  expect_equal(SPRDist(tr1, tr3, T), SPRDist(tr3, tr1, T))
  expect_equal(2L, length(SPRDist(trs12, tr3, T)))
  expect_equal(SPRDist(trs12, tr3, T), SPRDist(tr3, trs12, T))
  expect_equal(SPRDist(rev(trs123), trs123, T), SPRDist(trs123, rev(trs123), T))
  expect_equal(SPRDist(trs12, T))
  
})
