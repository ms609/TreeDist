context("Morphy: Tree search")

comb11 <- ape::read.tree(text="(a, (b, (c, (d, (e, (f, (g, (h, (i, (j, k))))))))));")
unrooted11 <- ape::read.tree(text="(a, b, (c, (d, (e, (f, (g, (h, (i, (j, k)))))))));")
data11 <- cbind(upper.tri(matrix(FALSE, 11, 11))[, 3:10], lower.tri(matrix(FALSE, 11, 11))[, 2:9])
rownames(data11) <- letters[1:11]
phy11 <- phangorn::phyDat(data11, type='USER', levels=c(FALSE, TRUE))
RootySwappers <- list(RootedTBRSwap, RootedSPRSwap, RootedNNISwap)

test_that("tree can be found", {
  set.seed(0)
  expect_error(TreeSearch(tree=unrooted11, dataset=phy11))
  expect_equal(comb11, TreeSearch(tree=RandomTree(phy11, 'a'), dataset=phy11,
               maxIter=2500, EdgeSwapper = RootedTBRSwap, verbosity=0))
  expect_equal(comb11, TreeSearch(RandomTree(phy11, 'a'), phy11, maxIter=2000, EdgeSwapper = RootedSPRSwap,
              verbosity=0))
  expect_equal(comb11, TreeSearch(RandomTree(phy11, 'a'), phy11, maxIter=2000, EdgeSwapper = RootedNNISwap, verbosity=0))
  expect_equal(comb11, Ratchet(RandomTree(phy11, 'a'), phy11, searchIter=300, searchHits = 20, swappers = RootySwappers, ratchHits=3, verbosity=0))
  expect_equal('multiPhylo', class(
    Ratchet(RandomTree(phy11, 'a'), phy11, searchIter=300, searchHits = 20,
    ratchHits=3, verbosity=0L, returnAll=TRUE)
    ))
  # expect_equal(Sectorial(RandomTree(phy11, 'a'), phy11, verbosity=-1), comb11) # TODO: Sectorial Search not working yet!
})


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

context("Implied weights: Tree search")
test_that("tree can be found", {
  set.seed(0)
  expect_error(IWTreeSearch(tree=unrooted11, dataset=phy11))
  expect_equal(comb11, IWTreeSearch(tree=RandomTree(phy11, 'a'), dataset=phy11,
                                  maxIter=2500, EdgeSwapper = RootedTBRSwap, verbosity=0))
  expect_equal(comb11, IWTreeSearch(RandomTree(phy11, 'a'), phy11, maxIter=2000, EdgeSwapper = RootedSPRSwap,
                                  verbosity=0))
  expect_equal(comb11, IWTreeSearch(RandomTree(phy11, 'a'), phy11, maxIter=2000, EdgeSwapper = RootedNNISwap, verbosity=0))
  expect_equal(comb11, IWRatchet(RandomTree(phy11, 'a'), phy11, searchIter=300, searchHits = 20, 
                                 swappers = RootySwappers, ratchHits=3, verbosity=0))
  expect_equal('multiPhylo', class(
    IWRatchet(RandomTree(phy11, 'a'), phy11, searchIter=300, searchHits = 20,
            ratchHits=3, verbosity=0L, returnAll=TRUE)
  ))
  # expect_equal(IWSectorial(RandomTree(phy11, 'a'), phy11, verbosity=-1), comb11) # TODO: Sectorial Search not working yet!
})

