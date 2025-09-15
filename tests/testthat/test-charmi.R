test_that("ExpSameCherries", {
  
  expect_simulated <- function(char, cherries) {
    exp <- ExpSameCherries(char, cherries)
    tt <- t.test(replicate(16000, sum(apply(matrix(sample(char, 2 * cherries, replace = FALSE), 2, cherries),
                       2, function(x) x[[1]] == x[[2]]))), mu = exp, conf.level = 0.997)
    conf <- tt[["conf.int"]]
    expect_gt(conf[[2]], exp)
    expect_lt(conf[[1]], exp)
  }
  
  char <- c(1,1,1,1,2,2,2)
  expect_simulated(char, 1)
  expect_simulated(char, 2)
  expect_simulated(char, 3)
  expect_error(ExpSameCherries(char, length(char)), "cherries can't be filled")
  
  char2 <- c(rep(1, 12), rep(2, 4), 4, rep(3, 11))
  expect_simulated(char2, 1)
  expect_simulated(char2, 3)
  expect_simulated(char2, 4)
  expect_simulated(char2, 6)
  expect_simulated(char2, 12)
  
  expect_emi <- function(char, tree) {
    
    nTip <- length(char)
    tree <- TreeTools::PectinateTree(nTip)
    nCherries <- Cherries(tree)
    treeH <- nTip * log2(nTip) - (2 * nCherries)
    expect_equal(treeH, CharH(tree))
    
    tab <- table(char, deparse.level = 0)
    charH <- (nTip * log2(nTip)) - sum(tab * log2(tab))
    expect_equal(charH, CharH(char))
    
    .ApproxEMI <- function(char, tree) {
      tree <- as.HPart(tree)
      char <- .MakeHPartMatch(char, tree)
      EMI_xptr(char, tree, 0.005, 36L) / log(2)
    }
    
    
    expect_equal(CharEMI(char, tree)[[1]], .ApproxEMI(char, tree)[[1]],
                 tolerance = 0.01)
    expect_equal(
      CharEMI(char, tree)[[1]],
      treeH + charH - (nTip * log2(nTip)) + 2 * ExpSameCherries(char, nCherries)
    )
  }
  
  
  expect_ami <- function(char, tree) {
    
    nTip <- length(char)
    tree <- TreeTools::PectinateTree(nTip)
    nCherries <- Cherries(tree)
    treeH <- nTip * log2(nTip) - (2 * nCherries)
    expect_equal(treeH, CharH(tree))
    
    tab <- table(char, deparse.level = 0)
    charH <- (nTip * log2(nTip)) - sum(tab * log2(tab))
    expect_equal(charH, CharH(char))
    
    .ApproxAMI <- function(char, tree) {
      tree <- as.HPart(tree)
      char <- .MakeHPartMatch(char, tree)
      AMI_xptr(char, tree, as.function(Mean), 0.005, 32L)
    }
    
    expect_equal(CharAMI(char, tree)[[1]], .ApproxAMI(char, tree)[[1]],
                 tolerance = 0.01)
    
    meanH <- charH
    jointH <- CharJH(char, tree)
    mi <- charH + treeH - jointH
    ejh <- CharEJH(char, tree)
    emi <- charH + treeH - ejh
    ami <- (mi - emi) / (meanH - emi)
    expect_equal(CharAMI(char, tree)[[1]], ami)
  }
  
  expect_emi(char2, TreeTools::PectinateTree(nTip))
  expect_emi(char2, TreeTools::BalancedTree(nTip))
  set.seed(1)
  expect_emi(char2, TreeTools::RandomTree(nTip))
})

test_that("CharMI works with simple trees", {
  Bal <- TreeTools::BalancedTree
  
  expect_equal(CharH(c(2, 2, 1, 1, 1)), 5 * Ntropy(2, 3))
  expect_equal(CharH(Bal(5)), 5 * log2(5) - 4)
  expect_equal(CharH(TreeTools::UnrootTree(Bal(5))), 5 * log2(5) - 4)
  expect_equal(CharH(TreeTools::PectinateTree(5)), 5 * log2(5) - 2)
  # !  - Fails: cherry at root not currently 'backstripped'
  expect_equal(CharH(TreeTools::UnrootTree(TreeTools::PectinateTree(5))), 5 * log2(5) - 4) 
  
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
  
  pec7 <- PectinateTree(7)
  expect_lt(CharMI(c(1,1,1,1,2,2,1), pec7), CharMI(c(1,1,1,1,1,2,2), pec7))
  
  # Mathematically correct test cases illustrate the shortcoming of the paradigm
  expect_equal(CharMI(c(1,1,1,2,2,1,1), pec7), CharMI(c(1,1,1,1,1,2,2), pec7))
  expect_gt(CharMI(c(1,2,2,1,1,1,1), pec7), CharMI(c(1,1,1,1,2,2,1), pec7))
  
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
  
  emi <- CharEMI(flatP, hp9, precision = 0.003)
  expect_lt(emi[[1]], h1)
  
  ami <- CharAMI(flatP, hp9, max, precision = 0.01)
  expect_equal(ami[[1]], (mi - emi[[1]]) / (max(h1, h2) - emi[[1]]),
               tolerance = 4 * 0.01)
  
  emAt <- attributes(emi)
  amAt <- attributes(ami)
  expect_equal(names(emAt), names(amAt))
  expect_equal(emAt[["ejh"]], amAt[["ejh"]], tolerance = emAt[["ejhSEM"]]
               + amAt[["ejhSEM"]])
  
  ami <- CharAMI(flatP, hp9)
  expect_equal((mi - emi[[1]]) / (h1 - emi[[1]]), 1)
  expect_equal(ami, structure(1, samples = 0, ejh = NA_real_, ejhVar = NA_real_,
                              ejhSD = NA_real_, ejhSEM = NA_real_,
                              sem = 0, precision = 0))
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
  
  # Check that the same seed returns the same output
  set.seed(1)
  samples2 <- replicate(256, CharAMI(ch, TreeTools::RandomTree(ch)))
  expect_equal(samples2, samples)
})
