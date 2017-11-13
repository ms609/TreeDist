library(ape)
library(testthat)

context("Morphy: tree search")
test_that("tree search finds shortest tree", {
  true_tree <- ape::read.tree(text = "(((((1,2),3),4),5),6);")
  malformed_tree <- ape::read.tree(text = "((((1,2),3),4),5,6);")
  dataset <- StringToPhyDat('110000 111000 111100', 1:6, byTaxon=FALSE)
  expect_error(TreeSearch(malformed_tree, dataset))
  start_tree <- RenumberTips(read.tree(text = "(((1, 6), 3), (2, (4, 5)));"), true_tree$tip.label)
  expect_equal(Fitch(start_tree, dataset), 6)
  morphyObj <- PhyDat2Morphy(dataset)
  on.exit(morphyObj <- UnloadMorphy(morphyObj))
  
  expect_equal(3, attr(TreeSearch(start_tree, dataset, EdgeSwapper=NNISwap, verbosity=0), 'score'),
               Fitch(true_tree, dataset))
  expect_equal(3, attr(TreeSearch(start_tree, dataset, EdgeSwapper=SPRSwap, verbosity=-1), 'score'),
               Fitch(true_tree, dataset))
  expect_equal(3, attr(TreeSearch(start_tree, dataset, EdgeSwapper=TBRSwap, verbosity=-1), 'score'),
               Fitch(true_tree, dataset))
  expect_equal(3, attr(TreeSearch(start_tree, dataset, EdgeSwapper=RootedNNISwap, verbosity=-1), 'score'),
               Fitch(true_tree, dataset))
  expect_equal(3, attr(TreeSearch(start_tree, dataset, EdgeSwapper=RootedSPRSwap, verbosity=-1), 'score'),
               Fitch(true_tree, dataset))
  expect_equal(3, attr(TreeSearch(start_tree, dataset, EdgeSwapper=RootedTBRSwap, verbosity=-1), 'score'),
               Fitch(true_tree, dataset))
  ratchetScore <- attr(Ratchet(start_tree, dataset, 
                  swappers=list(RootedTBRSwap, RootedSPRSwap, RootedNNISwap),
                  ratchIter=3, searchHits=5, verbosity=0), 'score')
  expect_equal(3, Fitch(true_tree, dataset), ratchetScore)
})
