test_that("CharMI works with real dataset", {
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
  #expect_equal(9 * Ntropy(table(ch)))
  expect_equal(CharH(bal9),
               NTip(bal9) * log2(NTip(bal9)) # Entropy of identifying each tip
               - 4 * (2 * log2(2)) # But we can't distinguish between cherries
               )
  expect_equal(CharH(tr), NTip(tr) * log2(NTip(tr)) - (3 * 2))
  expect_equal(CharH(StarTree(10)), 0)
  
  # Build HPart from tree, then relabel
  trPart <- as.HPart(tr)
  attr(trPart, "tip.label") <- seq_along(attr(trPart, "tip.label"))
  expect_equal(attr(chPart, "tip.label"), attr(trPart, "tip.label"))
  
  # Because of the difference in levels, this test should NOT pass (!)
  # expect_equal(HMI(chPart, trPart), SelfHMI(chPart))
  
  # Relabel tree first, then build HPart
  tree <- tr
  tree$tip.label <- seq_along(tree[["tip.label"]])
  treePart <- as.HPart(tree)
  treePart
  expect_equal(HMI(trPart, treePart), SelfHMI(treePart))
  expect_equal(HMI(chPart, trPart), HMI(chPart, treePart))
})

