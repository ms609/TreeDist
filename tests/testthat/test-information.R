library('TreeTools')

test_that("Entropy is calculated correctly", {
  expect_equal(1, Entropy(rep(0.5, 2)))
  expect_equal(2, Entropy(c(1/4, 1/4, 0, 1/4, 0, 1/4)))
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
  
  Test <- function(nTaxa, nContra) {
    nInSplit <- nTaxa / 2L
    split1 <- split2 <- c(rep(TRUE, nInSplit), rep(FALSE, nInSplit))
    flips <- nInSplit + ((1-nContra):nContra)
    nonFlips <- c(seq_len(nContra), nTaxa + 1 - seq_len(nContra))
    split2[flips] <- !split2[flips]
 
    expect_true(
      SharedPhylogeneticInfoSplits(as.Splits(split1[-flips]), as.Splits(split2[-flips]))
      >
      SharedPhylogeneticInfoSplits(as.Splits(split1[-nonFlips]), as.Splits(split2[-nonFlips]))
    )
  }
  
  Test(200, 2)
  Test(200, 1)
  Test(200, 10)
})


test_that("TreesConsistentWithTwoSplits works", {
  
  Test <- function(n, a, b, score) {
    logScore <- log(score)
    
    expect_equal(score, TreesConsistentWithTwoSplits(n, a, b))
    expect_equal(score, TreesConsistentWithTwoSplits(n, b, a))
    expect_equal(score, TreesConsistentWithTwoSplits(n, n - a, n - b))
    expect_equal(score, TreesConsistentWithTwoSplits(n, n - b, n - a))
    expect_equal(logScore, LnTreesConsistentWithTwoSplits(n, a, b))
    expect_equal(logScore, LnTreesConsistentWithTwoSplits(n, b, a))
    expect_equal(logScore, LnTreesConsistentWithTwoSplits(n, n - a, n - b))
    expect_equal(logScore, LnTreesConsistentWithTwoSplits(n, n - b, n - a))
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
  expect_equal(0.4199732, tolerance = 1e-6,
               MeilaMutualInformation(c(T,T,T,F,F), c(F,T,T,F,F)))
})

test_that("MeilaVariationOfInformation", {
  expect_equal(1L, MeilaVariationOfInformation(c(T,T,T,F,F,F), c(T,T,T,T,T,T)))
  expect_equal(0L, MeilaVariationOfInformation(c(T,T,T,F,F,F), c(T,T,T,F,F,F)))
  expect_equal(11.01955 / 6, MeilaVariationOfInformation(c(T,T,T,F,F,F), c(T,F,T,F,T,F)))
  expect_equal(7.219281 / 6, MeilaVariationOfInformation(c(F,T,T,T,T,T), c(T,T,T,T,T,F)))
})

test_that("SplitEntropy", {
  expect_equal(SplitEntropy(c(rep(TRUE, 5), rep(FALSE, 6)), c(rep(TRUE, 5), rep(FALSE, 6))),
               SplitEntropy(c(rep(TRUE, 5), rep(FALSE, 6))))
  expect_equal(c(H1 = 0.994, H2 = 0.994, H12 = 1.348, I = 0.6394, Hd = 0.709),
               SplitEntropy(c(rep(TRUE, 5), rep(FALSE, 6)), c(rep(TRUE, 6), rep(FALSE, 5))),
               tolerance = 0.001)
})
