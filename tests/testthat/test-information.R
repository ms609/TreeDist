context("Information.R")
library('TreeTools')

test_that("Entropy is calculated correctly", {
  expect_equal(1, Entropy(rep(0.5, 2)))
  expect_equal(2, Entropy(c(1/4, 1/4, 0, 1/4, 0, 1/4)))
})

test_that("Joint information calculated correctly", {
  # Identical splits: ABCDE:FGH, ABCDE:FGH
  expect_equal(-log2(315/10395), JointInformation(5, 0, 0, 3))
  expect_equal(-log2(315/10395), JointInformation(3, 0, 0, 5))
  expect_equal(-log2(315/10395), JointInformation(0, 5, 3, 0))
  expect_equal(-log2(315/10395), JointInformation(0, 3, 5, 0))
  
  # Agreeable splits: ABCDE:FGHI, ABC:DEFGHI
  expected_3204 <- SplitInformation(5, 4) + SplitInformation(3, 6) --log2(135/135135)
  expect_equal(expected_3204, JointInformation(3, 2, 0, 4))
  expect_equal(expected_3204, JointInformation(2, 3, 4, 0))
  expect_equal(expected_3204, JointInformation(0, 4, 3, 2))
  expect_equal(expected_3204, JointInformation(4, 0, 2, 3))
  
  # Perfect contradiction: AB:CDEFG, AC:BDEFG
  expect_equal(SplitInformation(2, 5) * 2, JointInformation(1, 1, 1, 4))
  
  # Compatible splits: AB:CDEFGH, CD:ABEFGH
  expected_0224 <- SplitInformation(2, 6) - SplitInformation(2, 5)
  expect_equal(expected_0224, JointInformation(0, 2, 2, 4))
  expect_equal(expected_0224, JointInformation(2, 0, 4, 2))
  expect_equal(expected_0224, JointInformation(2, 4, 0, 2))
  expect_equal(expected_0224, JointInformation(4, 2, 2, 0))
  
})

test_that("AllSplitPairings counted correctly", {
  expect_error(AllSplitPairings(3))
  for (n in 4:10) {
    totalSplits <- sum(choose(n, 2:(n-2)))
    expect_equal(totalSplits * totalSplits, sum(AllSplitPairings(n)))
  }
})

test_that("Removing contradictions improves scores", {
  # Imagine 200 taxa, split 
  #   ...00000111111.....
  #   ...00011001111.....
  # We can delete taxa 99 and 100 to produce matching splits of
  # 198=100:98, or we can delete taxa 99-102 to produce 196=96:96.
  
  Test <- function (nTaxa, nContra) {
    nInSplit <- nTaxa / 2L
    split1 <- split2 <- c(rep(TRUE, nInSplit), rep(FALSE, nInSplit))
    flips <- nInSplit + ((1-nContra):nContra)
    nonFlips <- c(seq_len(nContra), nTaxa + 1 - seq_len(nContra))
    split2[flips] <- !split2[flips]
    
    print(as.Splits(split1[-flips]), T)
    print(as.Splits(split2[-flips]), T)
    cat(MutualPhylogeneticInfoSplits(as.Splits(split1[-flips]), as.Splits(split2[-flips])))
    print(as.Splits(split1[-nonFlips]), T)
    print(as.Splits(split2[-nonFlips]), T)
    cat(MutualPhylogeneticInfoSplits(as.Splits(split1[-nonFlips]), as.Splits(split2[-nonFlips])))
    
    expect_true(
      MutualPhylogeneticInfoSplits(as.Splits(split1[-flips]), as.Splits(split2[-flips]))
      >
      MutualPhylogeneticInfoSplits(as.Splits(split1[-nonFlips]), as.Splits(split2[-nonFlips]))
    )
  }
  
  Test(10, 2)
  Test(10, 1)
  Test(10, 10)
  
  Test(200, 2)
  Test(200, 1)
  Test(200, 10)
})


test_that("TreesConsistentWithTwoSplits works", {
  
  Test <- function (n, a, b, score) {
    logScore <- log(score)
    
    expect_equal(score, TreesConsistentWithTwoSplits(n, a, b))
    expect_equal(score, TreesConsistentWithTwoSplits(n, b, a))
    expect_equal(score, TreesConsistentWithTwoSplits(n, n - a, n - b))
    expect_equal(score, TreesConsistentWithTwoSplits(n, n - b, n - a))
    expect_equal(logScore, LogTreesConsistentWithTwoSplits(n, a, b))
    expect_equal(logScore, LogTreesConsistentWithTwoSplits(n, b, a))
    expect_equal(logScore, LogTreesConsistentWithTwoSplits(n, n - a, n - b))
    expect_equal(logScore, LogTreesConsistentWithTwoSplits(n, n - b, n - a))
  }
  
  Test(8, 3, 0, 315)
  Test(8, 8, 0, NUnrooted(8))
  Test(8, 3, 3, TreesMatchingSplit(3, 5))
  Test(10, 5, 2, 1575)
  Test(9, 5, 3, 135)
  Test(8, 7, 3, 315)
})

test_that("MeilaMutualInformation", {
  expect_error(MeilaMutualInformation(c(T,T,T), c(T,T,T,T)))
  expect_equal(0, MeilaMutualInformation(c(T,T,T,F,F), c(F,F,F,F,F)))
})

test_that("MeilaVariationOfInformation", {
  expect_equal(6L, MeilaVariationOfInformation(c(T,T,T,F,F,F), c(T,T,T,T,T,T)))
  expect_equal(0, MeilaVariationOfInformation(c(T,T,T,F,F,F), c(T,T,T,F,F,F)))
  expect_equal(11.01955, MeilaVariationOfInformation(c(T,T,T,F,F,F), c(T,F,T,F,T,F)))
  expect_equal(7.219281, MeilaVariationOfInformation(c(F,T,T,T,T,T), c(T,T,T,T,T,F)))
})

test_that("SplitEntropy", {
  expect_equal(SplitEntropy(c(rep(TRUE, 5), rep(FALSE, 6)), c(rep(TRUE, 5), rep(FALSE, 6))),
               SplitEntropy(c(rep(TRUE, 5), rep(FALSE, 6))))
  expect_equal(c(h1=0.994, h2=0.994, jointH=1.348, i=0.6394, vI = 0.709),
               SplitEntropy(c(rep(TRUE, 5), rep(FALSE, 6)), c(rep(TRUE, 6), rep(FALSE, 5))),
               tolerance = 0.001)
})
