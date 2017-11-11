context("Profile Parsimony: Tree search")
test_that("Profile parsimony works in tree search", {
  sillyData <- lapply(1:22, function (i) c( rep(0, i - 1), rep(1, 22 - i), rep(1, 22 - i), rep(0, i - 1)))#, sample(2, 20, replace=TRUE)-1))
  names(sillyData) <- as.character(1:22)
  dataset <- PhyDat(sillyData, 0:1)
  readyData <- PrepareDataProfile(dataset, 12000)
  
  set.seed(0)
  rTree <- randomTree <- RandomTree(dataset, '1')
  expect_equal(FitchScore(rTree, readyData, TipsAreColumns), FitchScore(rTree, dataset, TipsAreNames))
  expect_equal(90, FitchScore(referenceTree, dataset, TipsAreNames))
  expect_equal(90, FitchScore(referenceTree, readyData, TipsAreColumns))
  
  expect_true(ProfileScore(rTree, readyData) > ProfileScore(referenceTree, readyData))

  quickTS <- TreeSearch(rTree, dataset, TreeScorer=FitchScore, EdgeSwapper=NNISwap, 
                        maxIter=10000, maxHits=40, verbosity=0)
  expect_equal(42, attr(quickTS, 'score'))
  quickFitch <- Ratchet(rTree, dataset, TreeScorer = FitchScore, suboptimal=max(PP_SUBOPTIMAL_VALUES),
                        ratchHits=3, searchHits=15, searchIter=500, ratchIter=500,
                        rearrangements='TBR', verbosity=-1)
  expect_equal(42, attr(quickFitch, 'score'))
                   
  quick <- Ratchet(rTree, readyData, TreeScorer = ProfileScore, returnAll = FALSE, rooted=TRUE,
                   ratchHits=5, searchHits=30, searchIter=100, ratchIter=50,
                   rearrangements='TBR', verbosity=-1)
  expect_equal(quick, quickFitch)
})
