context("Support.R")

test_that("Splits are counted correctly", {
  abCD <- matrix(c(FALSE, FALSE, TRUE, TRUE), dimnames=list(LETTERS[1:4], 2))
  ABcd <- matrix(c(TRUE, TRUE, FALSE, FALSE), dimnames=list(LETTERS[1:4], 1))
  AcBd <- matrix(c(TRUE, FALSE, TRUE, FALSE), dimnames=list(c('A', 'C', 'B', 'D'), 'X'))
  dABc <- matrix(c(FALSE, TRUE, TRUE, FALSE), dimnames=list(c('D', 'A', 'B', 'C'), NULL))
  
  expect_true(SplitsRepeated(abCD, ABcd))
  expect_true(SplitsRepeated(abCD, AcBd))
  expect_true(SplitsRepeated(abCD, dABc))
})

test_that("Jackknife supports are correct", {
  true_tree <-  ape::read.tree(text = "((((((A,B),C),D),E),F),out);")
  start_tree <- ape::read.tree(text = "(((((A,D),B),E),(C,F)),out);")
  dataset <- StringToPhyDat('1100000 1110000 1111000 1111100 1100000 1110000 1111000 1111100 1001000', 1:7, byTaxon=FALSE)
  names(dataset) <- c(LETTERS[1:6], 'out')
  set.seed(0)
  strict <- TreeSearch(start_tree, dataset)
  expect_equal(1, length(unique(list(true_tree), list(start_tree)))) # Right tree found
  jackTrees <- Jackknife(strict, dataset, resampleFreq=4/7, searchIter=200L, searchHits=7L, 
                         EdgeSwapper=RootedTBRSwap, jackIter=20L, verbosity=0L)
  # Note: one cause of failure could be a change in characters sampled, due to randomness
  expect_true(length(unique(jackTrees)) > 2L)
})

test_that("Node supports calculated correctly", {
  treeSample <- list(
    correct = ape::read.tree(text = "((((((A,B),C),D),E),F),out);"),
    swapFE  = ape::read.tree(text = "((((((A,B),C),D),F),E),out);"),
    DEClade = ape::read.tree(text = "(((((A,B),C),(D,E)),F),out);"),
    swapBC  = ape::read.tree(text = "((((((A,C),B),D),E),F),out);"),
    DbyA    = ape::read.tree(text = "((((((A,D),C),B),E),F),out);")
  )
  expect_equal(c('10'=4, '11'=4, '12'=4, '13'=3), 
               SplitFrequency(treeSample$correct, treeSample))
  
  # Internal nodes on each side of root
  balanced <- ape::read.tree(text="((D, (E, (F, out))), (C, (A, B)));")
  expect_equal(c('9'=4, '10'=4, '11'=4, '13'=3),
               SplitFrequency(balanced, treeSample))
  
})

test_that("Node support colours consistent", {
  expect_equal('red', SupportColour(NA))
  expect_equal('red', SupportColour(2))
  expect_equal('red', SupportColor(-2)) # Check alternative spelling 
  expect_equal('#ffffff00', SupportColour(1, show1=FALSE))
})
