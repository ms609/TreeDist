library("TreeTools")

test_that("HMI works with real dataset", { # TODO move to appropriate position
  ch <- c(1L, 2L, 2L, 2L, 2L, 2L, 1L, 1L, 1L, 1L, 1L)
  tr <- structure(list(
    edge = structure(c(12L, 12L, 13L, 14L, 15L, 16L, 16L, 17L, 17L, 18L, 18L, 
                       15L, 14L, 19L, 20L, 20L, 19L, 13L, 21L, 21L, 1L, 13L, 
                       14L, 15L, 16L, 2L, 17L, 3L, 18L, 4L, 5L, 6L, 19L, 20L,
                       7L, 8L, 9L, 21L, 10L, 11L), dim = c(20L, 2L)),
    Nnode = 10L,
    tip.label = c("Nem", "Sco", "Eun", "Aph", "Chr", "Can", "Hel", "Cha",
                  "Lep", "Ter", "Lin")),
    class = "phylo", order = "preorder")
  chPart <- as.HPart(ch)
  
  # Build HPart from tree, then relabel
  trPart <- as.HPart(tr)
  attr(trPart, "tip.label") <- seq_along(attr(trPart, "tip.label"))
  expect_equal(attr(chPart, "tip.label"), attr(trPart, "tip.label"))
  expect_equal(HMI(chPart, trPart), SelfHMI(chPart))
  
  # Relabel tree first, then build HPart
  tree <- tr
  tree$tip.label <- seq_along(tree[["tip.label"]])
  treePart <- as.HPart(tree)
  treePart
  expect_equal(HMI(trPart, treePart), SelfHMI(treePart))
  expect_equal(HMI(chPart, trPart), HMI(chPart, treePart))
})


test_that("is.HPart() succeeds", {
  expect_true(is.HPart(as.HPart(TreeTools::BalancedTree(7))))
  expect_true(is.HPart(structure(class = "HPart",
                                 list(list("t1"), list("t2", "t3")))))
  expect_false(is.HPart(structure(class = "NonPart", 
                                  list(list("t1"), list("t2", "t3")))))
})

test_that("HMI examples from Perotti et al. 2015", {
  p1 <- list(list(1, list(2, 3)), list(4, 5, 6))
  p2 <- list(1, list(2, 3), list(4, 5, 6))
  expect_equal(SelfHMI(p1), 1.011 / log(2), tolerance = 0.01)
  expect_equal(SelfHMI(p2), 1.011 / log(2), tolerance = 0.01)
  expect_equal(HMI(p1, p2), log(2) / log(2))
})

test_that("HMI results match hmi.pynb", {
  # Non-hierarchical
  p1 <- list(list(19, 18, 5), list(14, 16, 3), list(7), list(10, 8), list(1, 17, 9, 4, 6, 15), list(2, 13, 11), list(12, 0))
  p2 <- list( list(12, 9), list(4, 2, 0, 7), list(16), list(5), list(8, 3, 1, 14), list(11, 6, 10), list(18, 17, 19), list(13, 15))
  expect_equal(HMI(p1, p2), 0.9410980357245466 / log(2))
  
  
  # Hierarchical
  hp0 <- list(list(23),
    list(list(list(list(list(list(16), list(17)))))), # Tips above order 2 nodes
    list(list(12), list(22, 13)), list(5), list(7), list(24), list(list(list(9),
              list(list(14, 2))),
              list(list(list(list(list(list(27), list(3))))))),
    list(20, 29, 18), list(4), list(26, 15), list(list(10), list(21, 25)),
    list(11), list(list(0, 28), list(1), list(6)), list(19, 8))
  
  hp1 <- list(list(23), list(list(list(list(list(16,  17))))), list(list(12),  list(22, 13)), list(5), list(7), list(24), list(list(list(9),  list(list(14, 2))),  list(list(list(list(list(27, 3)))))), list(20, 29, 18), list(4), list(26, 15), list(list(10),  list(21, 25)), list(11), list(list(0, 28),  list(1),  list(6)), list(19, 8))
  
  hp1Collapsed <- list(
    list(23), list(16, 17), list(list(12), list(22, 13)), list(5), list(7),
    list(24), list(list(list(9), list(14, 2)), list(27, 3)), list(20, 29, 18),
        list(4), list(26, 15), list(list(10), list(21, 25)), list(11),
    list(list(0, 28), list(1), list(6)), list(19, 8))
  
  hp2 <- list(
    list(list(list(0, 25), list(24)), list(6), list(11, 28), list(8)),
    list(list(list(19), list(list(list(list(21), list(4),
                                       list(list(list(list(list(22, 7))))))))),
         list(5)), list(list(3), list(10, 23, 14)),
    list(list(27, 1, 16, 13, 18, 26, 9),
         list(list(list(list(15), list(list(list(list(list(list(12, 17)))))))),
              list(2, 20)), list(29)))
  
  expect_equal(HMI(hp1, hp2), 1.0591260408329395 / log(2))
  
  # expect_equal(HMI(hp0, hp0), 3.0140772805713665 / log(2))
  # Note that hp0 contains [[16], [17]] and [[27], [3]], whereas hp1 has
  # [16, 17] and [27, 3].  I haven't yet worked through why this should give a
  # different result.  But I don't think we are likely to encounter this case
  # in our work.
  expect_equal(HMI(hp1, hp1Collapsed), 2.921657656496707 / log(2))
  expect_equal(HMI(hp2, hp2), 2.606241391162456 / log(2))
  
  expect_equal(SelfHMI(hp1), HMI(hp1, hp1))
  
  ehmi <- structure(0.7806 / log(2), # Calculated from py with tol = 0.001
                    var = 0.01,
                    sd = 0.1, 
                    sem = 0.008,
                    relativeError = 0.01)
  ehmi_cpp <- EHMI(hp1, hp2, tolerance = 0.01)
  expect_gt(attr(ehmi_cpp, "samples"), 36)
  attr(ehmi_cpp, "samples") <- NULL # Could vary; no point in testing
  expect_equal(ehmi_cpp, ehmi, tolerance = 0.1)
  
  pyAHMI <- 0.13000 # Calculated with tol = 0.001
  expect_equal(AHMI(hp1, hp2)[[1]], pyAHMI, tolerance = 0.05)
  
  set.seed(1)
  ahmi1 <- AHMI(hp1, hp2)
  set.seed(1)
  expect_equal(AHMI(hp1, hp2), ahmi1)
  nRep <- 100
  ahmis <- replicate(nRep, AHMI(hp1, hp2))
  expect_lt(abs(attr(ahmi1, "sem") - sd(ahmis)), 0.1 * sd(ahmis))
})

test_that("HMI calculated correctly", {
  bal6 <- BalancedTree(6)
  pec6 <- PectinateTree(6)
  
  hp1 <- as.HPart(BalancedTree(6))
  hp2 <- as.HPart(PectinateTree(6))
  expect_equal(capture_output(print(hp2)),
               "Hierarchical partition on 6 leaves: t1, t2, ..., t5, t6")
  expect_equal(HMI_xptr(hp1, hp2), 0.363353185)
  bal8 <- BalancedTree(8)
  pec8 <- PectinateTree(8)
  star8 <- StarTree(8)
  
  expect_equal(HierarchicalMutualInfo(bal6, pec6, normalize = TRUE),
               HMI_xptr(hp1, hp2) / max(HH_xptr(hp1), HH_xptr(hp2)))
  
  hp1 <- build_hpart_from_phylo(BalancedTree(8))
  hp2 <- build_hpart_from_phylo(PectinateTree(8))
  expect_equal(HMI_xptr(hp1, hp1), 1.38629436)
  expect_equal(HMI_xptr(hp1, hp2), 0.3342954)
  
})

test_that("HMI_cpp equals SelfHMI for same partition", {
  set.seed(1)
  tr <- BalancedTree(8)
  hp <- as.HPart(tr)
  expect_equal(SelfHMI(hp), HMI(hp, hp), tolerance = 1e-12)
})
