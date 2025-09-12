test_that("CharMI works with simple trees", {
  Bal <- TreeTools::BalancedTree
  expect_equal(CharJH(c(2, 2, 1, 1, 1), c(1, 1, 2, 2, 2)),
               5 * Ntropy(2, 3))
  
  expect_equal(CharJH(c(2, 2, 1, 1, 1), c(1, 1, 1, 2, 2)),
               5 * Ntropy(2, 1, 2))
  
  expect_equal(CharJH(c(2, 2, 2, 1, 1), Bal(5)),
               CharH(Bal(5)))
  expect_equal(CharJH(c(2, 2, 1, 1, 1), Bal(5)),
               CharH(Bal(5)))
  
  bal9 <- Bal(9)
  expect_equal(CharH(bal9),
               NTip(bal9) * log2(NTip(bal9)) # Entropy of identifying each tip
               - 4 * (2 * log2(2)) # But we can't distinguish between cherries
  )
  expect_equal(CharH(TreeTools::StarTree(10)), 0)
  
  # Trees with >64 bits to test block logic
  
  expect_equal(CharJH(c(rep(2, 64), rep(1, 65)),
                      c(rep(1, 32), 2, rep(3, 32 + 64))),
               129 * Ntropy(32, 1, 31, 65))
  bal129 <- Bal(129)
  expect_equal(CharH(bal129),
               NTip(bal129) * log2(NTip(bal129)) # Entropy of identifying each tip
               - 128) # But we can't distinguish between cherries
})

test_that("CharAMI arithmetic checks out", {
  set.seed(1)
  
  flatP <- as.HPart(list(as.list(1:5), as.list(6:9)))
  hp9 <- as.HPart(TreeTools::BalancedTree(1:9))
  
  mi <- CharMI(flatP, hp9)
  expect_equal(mi, 0.99107606 * 9)
  
  h1 <- CharH(flatP)
  expect_equal(h1, 0.99107606 * 9)
  
  h2 <- CharH(hp9)
  expect_gt(h2, h1)
  
  emi <- CharEMI(flatP, hp9, precision = 0.003)[[1]]
  expect_lt(emi, h1)
  
  ami <- CharAMI(flatP, hp9, max, precision = 0.003)
  expect_equal(ami[[1]], (mi - emi) / (max(h1, h2) - emi), tolerance = 4 * 0.01)
  
  ami <- CharAMI(flatP, hp9)
  expect_equal((mi - emi) / (h1 - emi), 1)
  expect_equal(ami, structure(1, sem = 0))
})

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
  expect_equal(CharH(ch), length(ch) * Ntropy(table(ch)))
  
  nTip <- TreeTools::NTip(tr)
  expect_equal(CharH(tr), nTip * log2(nTip) - (3 * 2)) # Three cherries.
  
  expect_lt(CharJH(ch, tr), CharH(ch) + CharH(tr))
  expect_gte(CharJH(ch, tr), CharH(tr))
  
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

test_that("AHMI succeeds with CharMI", {
  ch <- c(1L, 2L, 2L, 2L, 2L, 2L, 1L, 1L, 1L, 1L, 1L)
  tr <- TreeTools::BalancedTree(11)
  XLogX <- function(x) x * log2(x)
  expect_equal(CharH(ch), XLogX(11) - XLogX(6) - XLogX(5))
  expect_equal(CharH(tr), XLogX(11) - 2 * 4)
  
  expect_error(expect_warning(CharJH(tr, ch), "Char.* a tree"), "t11 missing")
  expect_error(expect_warning(CharMI(tr, ch), "Char.* a tree"), "t11 missing")
  expect_error(expect_warning(CharEJH(tr, ch), "Char.* a tree"), "t11 missing")
  expect_error(expect_warning(CharEMI(tr, ch), "Char.* a tree"), "t11 missing")
  expect_error(expect_warning(CharAMI(tr, ch), "Char.* a tree"), "t11 missing")
  
  set.seed(1)
  expect_equal(CharH(ch) + CharH(tr) - CharJH(ch, tr), CharMI(ch, tr))
  expect_equal(CharH(ch) + CharH(tr) - CharEJH(ch, tr, precision = 0.001)[[1]],
               CharEMI(ch, tr, precision = 0.001)[[1]], tolerance = 0.02)
  
  expect_lt(CharAMI(ch, tr, max), CharAMI(ch, tr, min))
  expect_lt(CharAMI(c(1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1), tr, min)[[1]], 0)
  expect_equal(CharAMI(c(1, 1, rep(0, 9)), tr)[[1]], 1)
  expect_equal(CharAMI(c(rep(TRUE, 6), rep(FALSE, 5)), tr, min)[[1]], 1)
})

test_that("AHMI returns zero for random trees", {
  ch <- c(1L, 2L, 2L, 2L, 2L, 2L, 1L, 1L, 1L, 1L, 1L)
  tree <- TreeTools::NJTree(TreeTools::StringToPhyDat("12222211111"))
  expect_gt(CharAMI(ch, tree), 0)
  set.seed(1)
  samples <- replicate(256, CharAMI(ch, TreeTools::RandomTree(ch)))
  ci <- t.test(samples, mu = 0, conf.level = 0.997)$conf.int
  expect_lt(ci[[1]], 0)
  expect_gt(ci[[2]], 0)
})
