library(ape)
library(testthat)

context("Morphy: tree search")
test_that("tree search finds shortest tree", {
  true_tree <- ape::read.tree(text = "(((((1,2),3),4),5),6);")
  malformed_tree <- ape::read.tree(text = "((((1,2),3),4),5,6);")
  dataset <- StringToPhyDat('110000 111000 111100', 1:6, byTaxon=FALSE)
  expect_error(TreeSearch(malformed_tree, dataset))
  start_tree <- RenumberTips(read.tree(text = "(((1, 6), 3), (2, (4, 5)));"), true_tree$tip.label)
  expect_equal(InapplicableFitch(start_tree, dataset), 6)
  morphyObj <- PhyDat2Morphy(dataset)
  on.exit(morphyObj <- UnloadMorphy(morphyObj))
  
 ## TODO FIX VERBOISTY
  expect_equal(3, attr(MorphySearch(start_tree, dataset, Rearrange=NNICore, verbosity=-1), 'score'),
               InapplicableFitch(true_tree, dataset))
  expect_equal(3, attr(MorphySearch(start_tree, dataset, Rearrange=SPRCore, verbosity=-1), 'score'),
               InapplicableFitch(true_tree, dataset))
  expect_equal(3, attr(MorphySearch(start_tree, dataset, Rearrange=TBRCore, verbosity=-1), 'score'),
               InapplicableFitch(true_tree, dataset))
  expect_equal(3, attr(MorphySearch(start_tree, dataset, Rearrange=RootedNNICore, verbosity=-1), 'score'),
               InapplicableFitch(true_tree, dataset))
  expect_equal(3, attr(MorphySearch(start_tree, dataset, Rearrange=RootedSPRCore, verbosity=-1), 'score'),
               InapplicableFitch(true_tree, dataset))
  expect_equal(3, attr(MorphySearch(start_tree, dataset, Rearrange=RootedTBRCore, verbosity=-1), 'score'),
               InapplicableFitch(true_tree, dataset))
  ratchetScore <- attr(MorphyRatchet(start_tree, dataset, 
                  rearrangements=list(RootedTBRCore, RootedSPRCore, RootedNNICore),
                  k=3, maxHits=5, verbosity=0), 'score')
  expect_equal(3, InapplicableFitch(true_tree, dataset), ratchetScore)
})
