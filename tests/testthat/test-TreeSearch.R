context("TreeSearch.R")

comb11 <- ape::read.tree(text="(a, (b, (c, (d, (e, (f, (g, (h, (i, (j, k))))))))));")
unrooted11 <- ape::read.tree(text="(a, b, (c, (d, (e, (f, (g, (h, (i, (j, k)))))))));")
data11 <- cbind(upper.tri(matrix(FALSE, 11, 11))[, 3:10], lower.tri(matrix(FALSE, 11, 11))[, 2:9])
rownames(data11) <- letters[1:11]
phy11 <- phangorn::phyDat(data11, type='USER', levels=c(FALSE, TRUE))
RootySwappers <- list(RootedTBRSwap, RootedSPRSwap, RootedNNISwap)

test_that("tree can be found", {
  suppressWarnings(RNGversion("3.5.0")) # Until we can require R3.6.0
  set.seed(0)
  random11 <- RandomTree(phy11, 'a')
  expect_error(TreeSearch(unrooted11, dataset=phy11))
  expect_equal(TreeSearch(random11, dataset=phy11, maxIter=250, 
                          EdgeSwapper=RootedTBRSwap, verbosity=0L), comb11)
  expect_equal(TreeSearch(random11, dataset=phy11, maxIter=250, EdgeSwapper=AllTBR, 
                          stopAtPeak=TRUE, stopAtPlateau=10L, verbosity=0L), comb11)
  expect_equal(TreeSearch(random11, phy11, maxIter=400,
                                  EdgeSwapper=RootedSPRSwap, verbosity=0L), comb11)
  expect_equal(TreeSearch(RandomTree(phy11, 'a'), phy11, maxIter=200, EdgeSwapper = RootedNNISwap, verbosity=0), comb11)
  expect_equal(comb11, Ratchet(random11, phy11, searchIter=10, searchHits = 5,
                               swappers = RootySwappers, ratchHits=3, verbosity=0))
#  expect_equal(SectorialSearch(RandomTree(phy11, 'a'), phy11, verbosity=-1), comb11) # TODO: Sectorial Search not working yet!
})

test_that("tree search finds shortest tree", {
  true_tree <- ape::read.tree(text = "(((((1,2),3),4),5),6);")
  malformed_tree <- ape::read.tree(text = "((((1,2),3),4),5,6);")
  dataset <- StringToPhyDat('110000 111000 111100', 1:6, byTaxon=FALSE)
  expect_error(TreeSearch(malformed_tree, dataset))
  start_tree <- RenumberTips(ape::read.tree(text = "(((1, 6), 3), (2, (4, 5)));"), true_tree$tip.label)
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

test_that("Implied weights: Tree search", {
  suppressWarnings(RNGversion("3.5.0")) # Until we can require R3.6.0
  set.seed(0)
  expect_error(IWTreeSearch(tree=unrooted11, dataset=phy11))
  expect_equal(comb11, IWTreeSearch(RandomTree(phy11, 'a'), phy11, maxIter=100, 
                                    EdgeSwapper = RootedTBRSwap, verbosity=0))
  
  expect_equal(comb11, IWTreeSearch(RandomTree(phy11, 'a'), phy11, maxIter=200,
                                    EdgeSwapper = RootedSPRSwap, verbosity=0))
  
  expect_equal(comb11, IWTreeSearch(TBR(TBR(TBR((comb11)))), phy11, maxIter=100, 
                                    EdgeSwapper = RootedNNISwap, verbosity=0))
  
  expect_equal(comb11, IWRatchet(RandomTree(phy11, 'a'), phy11, searchIter=8,
                                 searchHits = 3, swappers = RootySwappers, 
                                 ratchHits = 3, verbosity=0))
  
  expect_equal('multiPhylo', class(
    IWRatchet(tree=RandomTree(phy11, 'a'), dataset=phy11, concavity=4,
              searchIter = 5, searchHits = 2,
              ratchHits = 2, verbosity=0L, returnAll=TRUE)
  ))
  # expect_equal(IWSectorial(RandomTree(phy11, 'a'), phy11, verbosity=-1), comb11) # TODO: Sectorial Search not working yet!
})


test_that("Profile parsimony works in tree search", {
  sillyData <- lapply(1:22, function (i) c( rep(0, i - 1), rep(1, 22 - i), rep(1, 22 - i), rep(0, i - 1)))#, sample(2, 20, replace=TRUE)-1))
  names(sillyData) <- as.character(1:22)
  dataset <- PhyDat(sillyData)
  readyData <- PrepareDataProfile(dataset, 12000, warn=FALSE)
  
  suppressWarnings(RNGversion("3.5.0")) # Until we can require R3.6.0
  set.seed(0)
  
  rTree <- randomTree <- RandomTree(dataset, '1')
  expect_equal(Fitch(rTree, readyData), Fitch(rTree, dataset))
  expect_equal(90, Fitch(referenceTree, dataset), Fitch(referenceTree, readyData))
  expect_true(ProfileScore(rTree, readyData) > ProfileScore(referenceTree, readyData))
  
  quickTS <- TreeSearch(rTree, dataset, TreeScorer=MorphyLength, EdgeSwapper=RootedNNISwap, 
                        maxIter=1600, maxHits=40, verbosity=0)
  expect_equal(42L, attr(quickTS, 'score'))
  
  quickFitch <- Ratchet(rTree, dataset, TreeScorer=MorphyLength, suboptimal=2, 
                        swappers=RootySwappers, ratchHits=3, searchHits=15, 
                        searchIter=100, ratchIter=500,
                        verbosity=0L)
  expect_equal(42, attr(quickFitch, 'score'))
  
  quickProf <- ProfileRatchet(rTree, readyData, 
                              swappers = RootySwappers,
                              BootstrapSwapper = RootedSPRSwap,
                              ratchIter=30L, ratchHits=3L, searchIter=30L, searchHits=3L,
                              verbosity=0L)
  expect_equal(quickProf, quickFitch)
})
