context('tree_distance.R')

test_that("Split combatibility is correctly established", {
  expect_true(SplitsCompatible(as.logical (c(0,0,1,1,0)), 
                               as.logical(c(0,0,1,1,0))))
  expect_true(SplitsCompatible(as.logical (c(0,0,1,1,0)), 
                               !as.logical(c(0,0,1,1,0))))
  expect_true(SplitsCompatible(as.logical (c(0,0,1,1,0)), 
                               as.logical(c(1,0,1,1,0))))
  expect_true(SplitsCompatible(!as.logical(c(0,0,1,1,0)), 
                               as.logical(c(1,0,1,1,0))))
  expect_false(SplitsCompatible(as.logical(c(0,0,1,1,0)), 
                                as.logical(c(1,1,0,1,0))))
})

methodsToTest <- list(
  MutualPhylogeneticInfo,
  VariationOfPhylogeneticInfo,
  MutualMatchingSplitInfo,
  VariationOfMatchingSplitInfo,
  MutualClusteringInfo,
  VariationOfClusteringInfo,
  NyeTreeSimilarity,
  JaccardRobinsonFoulds,
  MatchingSplitDistance,
  RobinsonFoulds,
  RobinsonFouldsInfo,
  KendallColijn
)

# Labels in different order to confound Tree2Splits
treeSym8 <- ape::read.tree(text='((e, (f, (g, h))), (((a, b), c), d));')
treeBal8 <- ape::read.tree(text='(((e, f), (g, h)), ((a, b), (c, d)));')
treeOpp8 <- ape::read.tree(text='(((a, f), (c, h)), ((g, b), (e, d)));')

treeCat8 <- ape::read.tree(text='((((h, g), f), e), (d, (c, (b, a))));')
treeTac8 <- ape::read.tree(text='((((e, c), g), a), (h, (b, (d, f))));')

treeAb.Cdefgh <- ape::read.tree(text='((a, b), (c, d, e, f, g, h));')
treeAbc.Defgh <- ape::read.tree(text='((a, b, c), (d, e, f, g, h));')
treeAcd.Befgh <- ape::read.tree(text='((a, c, d), (b, e, f, g, h));')
treeAbcd.Efgh <- ape::read.tree(text='((a, b, c, d), (e, f, g, h));')
treeTwoSplits <- ape::read.tree(text="(((a, b), c, d), (e, f, g, h));")

test_that('Bad labels cause error', {
  treeBadLabel8 <- ape::read.tree(text='((a, b, c, D), (e, f, g, h));')
  lapply(methodsToTest, function(Func) 
    expect_error(Func(treeSym8, treeBadLabel8)))
})

test_that('Size mismatch causes error', {
  treeSym7 <- ape::read.tree(text='((e, (f, g)), (((a, b), c), d));')
  lapply(methodsToTest, function(Func) 
    expect_error(Func(treeSym8, treeSym7)))
})

test_that('Metrics handle polytomies', {
  polytomy8 <- ape::read.tree(text='(a, b, c, d, e, f, g, h);')
  lapply(list(MutualPhylogeneticInfo, MutualClusteringInfo,
              MatchingSplitDistance, NyeTreeSimilarity),
         function (Func) expect_equal(0, Func(treeSym8, polytomy8)))
})

test_that('Output dimensions are correct', {
  Test <- function (Func) {
    answer <- 
    matrix(c(Func(treeSym8, treeSym8),      Func(treeBal8, treeSym8),
             Func(treeSym8, treeAbc.Defgh), Func(treeBal8, treeAbc.Defgh),
             Func(treeSym8, treeAbcd.Efgh), Func(treeBal8, treeAbcd.Efgh)),
           2L, 3L, dimnames=list(c('sym', 'bal'), c('sym', 'abc', 'abcd')))
    expect_equal(answer, Func(list(sym=treeSym8, bal=treeBal8), 
                              list(sym=treeSym8, abc=treeAbc.Defgh,
                                   abcd=treeAbcd.Efgh)))
  }
  lapply(methodsToTest, Test)
})

test_that('Mutual Phylogenetic Info is correctly calculated', {
  expect_equal(22.53747, tolerance=1e-05,
               MutualPhylogeneticInfo(treeSym8, treeSym8, normalize = FALSE))
  expect_equal(1, tolerance = 1e-05,
               MutualPhylogeneticInfo(treeSym8, treeSym8, normalize = TRUE))
  expect_equal(13.75284, MutualPhylogeneticInfo(treeSym8, treeBal8), tolerance=1e-05)
  expect_equal(VariationOfPhylogeneticInfo(treeSym8, treeAcd.Befgh),
               VariationOfPhylogeneticInfo(treeAcd.Befgh, treeSym8), tolerance=1e-05)
  expect_equal(0, VariationOfPhylogeneticInfo(treeSym8, treeSym8, normalize = TRUE))
  infoSymBal <- PartitionInfo(treeSym8) + PartitionInfo(treeBal8)
  expect_equal(infoSymBal - 13.75284 - 13.75284, tolerance = 1e-05,
    VariationOfPhylogeneticInfo(treeSym8, treeBal8, normalize = TRUE) * infoSymBal)
  expect_equal(22.53747 + MutualPhylogeneticInfo(treeAcd.Befgh, treeAcd.Befgh) - 
                 (2 * MutualPhylogeneticInfo(treeSym8, treeAcd.Befgh)), 
               VariationOfPhylogeneticInfo(treeSym8, treeAcd.Befgh), tolerance=1e-05)
  expect_equal(-log2(945/10395), MutualPhylogeneticInfo(treeSym8, treeAb.Cdefgh))
  expect_equal(22.53747 + MutualPhylogeneticInfo(treeBal8, treeBal8) - 13.75284 - 13.75284, 
               VariationOfPhylogeneticInfo(treeSym8, treeBal8), tolerance=1e-05)
  expect_equal(-log2(945/10395), MutualPhylogeneticInfo(treeSym8, treeAb.Cdefgh))
  expect_equal(-log2(315/10395), MutualPhylogeneticInfo(treeSym8, treeAbc.Defgh))
  expect_equal(0, VariationOfPhylogeneticInfo(treeSym8, treeSym8))
  expect_equal(PartitionInfo(treeSym8) - PartitionInfo(treeAcd.Befgh),
               VariationOfPhylogeneticInfo(treeSym8, treeAbc.Defgh))
  
  
  # Test symmetry of small vs large splits
  expect_equal(MutualPhylogeneticInfo(treeSym8, treeAbc.Defgh),
               MutualPhylogeneticInfo(treeAbc.Defgh, treeSym8))
  expect_equal(-log2(225/10395), MutualPhylogeneticInfo(treeSym8, treeAbcd.Efgh))
  expect_equal(-log2(225/10395) - log2(945/10395),
               MutualPhylogeneticInfo(treeSym8, treeTwoSplits))
  expect_equal(SplitMutualInformation(8, 4, 3),
               MutualPhylogeneticInfo(treeTwoSplits, treeAbc.Defgh))
  expect_equal(SplitInformation(4, 4) + SplitInformation (3, 5) - 
               (2 * SplitMutualInformation(8, 4, 3)),
               SplitVariationOfInformation(8, 4, 3), tolerance=1e-07)
  
  expect_equal(MutualPhylogeneticInfo(treeSym8, list(treeSym8, treeBal8)), 
               MutualPhylogeneticInfo(list(treeSym8, treeBal8), treeSym8))
})

test_that('MutualMatchingSplitInfo is correctly calculated', {
  BinaryToSplit <- function (binary) matrix(as.logical(binary))
  expect_equal(MutualMatchingSplitInfoSplits(
    BinaryToSplit(c(1, 1, 0, 0, 0, 0, 0, 0)),
    BinaryToSplit(c(0, 0, 1, 1, 0, 0, 0, 0))
    ), MutualMatchingSplitInfoSplits(
    BinaryToSplit(c(0, 0, 0, 0, 0, 0, 1, 1)),
    BinaryToSplit(c(0, 0, 1, 1, 0, 0, 0, 0))
    ))
  
  MutualMatchingSplitInfoSplits(BinaryToSplit(c(1, 1, 1, 1, 0, 0, 0, 0)),
                         BinaryToSplit(c(1, 0, 1, 0, 1, 0, 1, 0)))
  expect_equal(MutualPhylogeneticInfo(treeSym8, treeSym8),
               MutualMatchingSplitInfo(treeSym8, treeSym8), tolerance=1e-05)
  expect_equal(MutualMatchingSplitInfo(treeAb.Cdefgh, treeAbc.Defgh),
               MutualMatchingSplitInfo(treeAbc.Defgh, treeAb.Cdefgh))
  expect_equal(MutualMatchingSplitInfo(treeAbcd.Efgh, treeAb.Cdefgh),
               MutualMatchingSplitInfo(treeAb.Cdefgh, treeAbcd.Efgh))
  expect_equal(-(TreeSearch::LogTreesMatchingSplit(2, 5) - LnUnrooted.int(7)) / 
                 log(2), 
               MutualMatchingSplitInfo(treeAb.Cdefgh, treeAbc.Defgh))
  expect_true(MutualMatchingSplitInfo(treeSym8, treeBal8) > 
                MutualMatchingSplitInfo(treeSym8, treeOpp8))
  expect_equal(0, VariationOfMatchingSplitInfo(treeSym8, treeSym8))
})

test_that("Mutual Phylogenetic Information is correctly estimated", {
  exp <- ExpectedVariation(treeSym8, treeAbc.Defgh, samples=100)
  tol <- exp[, 'Std. Err.'] * 2
  # Expected values calculated with 10k samples
  expect_equal(exp['MutualPhylogeneticInfo', 'Estimate'], 
               1.166, tolerance=tol[1])
  expect_equal(exp['MutualMatchingSplitInfo', 'Estimate'], 
               3.096, tolerance=tol[2])
  expect_equal(exp['VariationOfPhylogeneticInfo', 'Estimate'], 
               25.250, tolerance=tol[3])
  expect_equal(exp['VariationOfMatchingSplitInfo', 'Estimate'], 
               21.388, tolerance=tol[4])
  expect_equal(exp[, 'sd'], exp[, 'Std. Err.'] * sqrt(exp[, 'n']))
})

test_that('Clustering information is correctly calculated', {
  expect_equal(ClusteringInfo(treeSym8), MutualClusteringInfo(treeSym8, treeSym8),
               tolerance=1e-05)
  expect_equal(1, MutualClusteringInfo(treeSym8, treeSym8, normalize = TRUE))
  expect_true(MutualClusteringInfo(treeSym8, treeBal8, normalize = pmin) >
                MutualClusteringInfo(treeSym8, treeBal8, normalize = pmax))
  expect_equal(ClusteringInfo(treeSym8) + ClusteringInfo(treeBal8) -
                 (2 * MutualClusteringInfo(treeBal8, treeSym8))
               , VariationOfClusteringInfo(treeSym8, treeBal8), tolerance=1e-05)
  expect_equal(MutualClusteringInfo(treeAb.Cdefgh, treeAbc.Defgh),
               MutualClusteringInfo(treeAbc.Defgh, treeAb.Cdefgh),
               tolerance=1e-05)
})

test_that('Matching Split Distance is correctly calculated', {
  expect_equal(0L, MatchingSplitDistance(treeSym8, treeSym8))
  expect_equal(1L, MatchingSplitDistance(treeAb.Cdefgh, treeAbc.Defgh))
  expect_equal(2L, MatchingSplitDistance(treeAb.Cdefgh, treeAbcd.Efgh))
  
  # Invariant to tree description order
  sq_pectinate <- ape::read.tree(text='((((((1, 2), 3), 4), 5), 6), (7, (8, (9, (10, 11)))));')
  shuffle1 <- ape::read.tree(text='(((((1, 5), 2), 6), (3, 4)), ((8, (7, 9)), (10, 11)));')
  shuffle2 <- ape::read.tree(text='(((8, (7, 9)), (10, 11)), ((((1, 5), 2), 6), (3, 4)));')
  expect_equal(MatchingSplitDistance(shuffle1, sq_pectinate),
               MatchingSplitDistance(sq_pectinate, shuffle1))
  expect_equal(0L, MatchingSplitDistance(shuffle1, shuffle2))
  expect_equal(MatchingSplitDistance(shuffle1, sq_pectinate),
               MatchingSplitDistance(shuffle2, sq_pectinate))
})

test_that('NyeTreeSimilarity is correctly calculated', {
  expect_equal(5L, NyeTreeSimilarity(treeSym8, treeSym8))
  expect_equal(c(3.8, 5), NyeTreeSimilarity(treeSym8, list(treeBal8, treeSym8)))
  expect_equal(2 / 3, NyeTreeSimilarity(treeAb.Cdefgh, treeAbc.Defgh))
  expect_equal(2 * (1 / 3), NyeTreeSimilarity(treeAb.Cdefgh, treeAbc.Defgh,
                                        similarity = FALSE))
  expect_equal(1, NyeTreeSimilarity(treeSym8, treeSym8, normalize = TRUE))
  #TODO: Validate expected value
  expect_equal(1L, NyeTreeSimilarity(treeSym8, treeAbcd.Efgh, 
                                     normalize = FALSE))
  expect_equal(1L / 5L, NyeTreeSimilarity(treeSym8, treeAbcd.Efgh, 
                                          normalize = TRUE))
  expect_true(NyeTreeSimilarity(treeSym8, treeBal8) > 
                NyeTreeSimilarity(treeSym8, treeOpp8))
})


test_that('Jaccard RF extremes tend to equivalent functions', {
  expect_equal(JaccardRobinsonFoulds(treeSym8, list(treeBal8, treeSym8),
                                     similarity = TRUE, k = 1L, arboreal = FALSE),
               NyeTreeSimilarity(treeSym8, list(treeBal8, treeSym8)) * 2L)
  
  expect_equal(JaccardRobinsonFoulds(treeSym8, list(treeBal8, treeSym8),
                                     similarity=FALSE, k=Inf),
               RobinsonFoulds(treeSym8, list(treeBal8, treeSym8)))
})

test_that('Jaccard RF is correctly calculated', {
  expect_equal(5L * 2L, JaccardRobinsonFoulds(treeSym8, treeSym8,
                                         k = 2, similarity = TRUE))
  expect_equal(c(3.32, 5) * 2L, 
               JaccardRobinsonFoulds(treeSym8, list(treeBal8, treeSym8),
                                     similarity = TRUE, k = 2))
  expect_equal(2 * 2, 3 * JaccardRobinsonFoulds(treeAb.Cdefgh, treeAbc.Defgh,
                                            similarity = TRUE))
  expect_equal(1, JaccardRobinsonFoulds(treeSym8, treeSym8,
                                        similarity = TRUE, normalize = TRUE))
  expect_equal(0, JaccardRobinsonFoulds(treeSym8, treeSym8,
                                        similarity = FALSE, normalize = TRUE))
  #TODO: Validate expected value
  expect_equal(1L * 2L, 
               JaccardRobinsonFoulds(treeSym8, treeAbcd.Efgh, similarity = TRUE,
                                     normalize = FALSE, k = 2))
  expect_equal(1L * 2L / 6L, 
               JaccardRobinsonFoulds(treeSym8, treeAbcd.Efgh, similarity = TRUE,
                                     normalize = TRUE, k = 2))
  expect_true(JaccardRobinsonFoulds(treeSym8, treeBal8, k = 2) < 
                JaccardRobinsonFoulds(treeSym8, treeOpp8, k = 2))
  expect_true(JaccardRobinsonFoulds(treeSym8, treeBal8, k = 3L) <
                JaccardRobinsonFoulds(treeSym8, treeBal8, k = 4L))
  expect_true(JaccardRobinsonFoulds(treeCat8, treeTac8, arboreal = FALSE) <
              JaccardRobinsonFoulds(treeCat8, treeTac8, arboreal = TRUE))
})

test_that('RobinsonFoulds is correctly calculated', {
  RF <- function (tree1, tree2) {
    suppressMessages(phangorn::RF.dist(tree1, tree2))
  }
  RFTest <- function (tree1, tree2) {
    expect_equal(RF(tree1, tree2), RobinsonFoulds(tree1, tree2))
  }
  RFTest(treeSym8, treeSym8)
  RFTest(treeBal8, treeSym8)
  expect_equal(c(4, 0), RobinsonFoulds(treeSym8, list(treeBal8, treeSym8)))
  RFTest(treeAb.Cdefgh, treeAbc.Defgh)
  expect_equal(0, RobinsonFoulds(treeSym8, treeSym8, normalize = TRUE))
  expect_equal(4L / 6L, 
               RobinsonFoulds(treeSym8, treeAbcd.Efgh, normalize = TRUE))
  RFTest(treeSym8, treeOpp8)
})


test_that('Robinson Foulds Info is correctly calculated', {
  expect_equal(22.53747 * 2L, tolerance=1e-05,
               RobinsonFouldsInfo(treeSym8, treeSym8, similarity = TRUE,
                                  normalize = FALSE))
  expect_equal(0, tolerance = 1e-05,
               RobinsonFouldsInfo(treeSym8, treeSym8, normalize = TRUE))
  expect_equal(1, tolerance = 1e-05,
               RobinsonFouldsInfo(treeSym8, treeSym8, similarity = TRUE, 
                                  normalize = TRUE))
  expect_equal(24.9, tolerance=0.01, 
               RobinsonFouldsInfo(treeSym8, treeBal8, similarity = TRUE))
  expect_equal(PartitionInfo(treeSym8) + PartitionInfo(treeBal8) -
                 RobinsonFouldsInfo(treeSym8, treeBal8, similarity = FALSE),
               RobinsonFouldsInfo(treeSym8, treeBal8, similarity = TRUE))
  expect_equal(-log2(945/10395) * 2,
               RobinsonFouldsInfo(treeSym8, treeAb.Cdefgh, similarity = TRUE))
  expect_equal(-log2(945/10395) * 2, 
               RobinsonFouldsInfo(treeSym8, treeAb.Cdefgh, similarity = TRUE))
  expect_equal(-log2(315/10395) * 2, 
               RobinsonFouldsInfo(treeSym8, treeAbc.Defgh, similarity = TRUE))
  
  # Test symmetry of small vs large splits
  expect_equal(RobinsonFouldsInfo(treeSym8, treeAbc.Defgh),
               RobinsonFouldsInfo(treeAbc.Defgh, treeSym8))
  expect_equal(-log2(225/10395) * 2, 
               RobinsonFouldsInfo(treeSym8, treeAbcd.Efgh, similarity = TRUE))
  expect_equal((-log2(225/10395) - log2(945/10395)) * 2,
               RobinsonFouldsInfo(treeSym8, treeTwoSplits, similarity = TRUE))
  expect_equal(RobinsonFouldsInfo(treeSym8, list(treeSym8, treeBal8)), 
               RobinsonFouldsInfo(list(treeSym8, treeBal8), treeSym8))
})


test_that('Kendall-Colijn distance is correctly calculated', {
  
    # Expected values calculated using treespace::treeDist(treeSym8, treeBal8)
  expect_equal(2.828427, KendallColijn(treeSym8, treeBal8), tolerance=1e-06)
  expect_equal(2.828427, KendallColijn(treeCat8, treeBal8), tolerance=1e-06)
  expect_equal(7.211103, KendallColijn(treeSym8, treeOpp8), tolerance=1e-06)
  expect_equal(matrix(c(0L, 8L), nrow=2, ncol=2, byrow=TRUE), 
               KendallColijn(list(treeSym8, treeCat8), list(treeCat8, treeTac8)), tolerance=1e-06)
  expect_equal(8L, KendallColijn(treeCat8, treeTac8), tolerance=1e-06)
  expect_equal(0L, KendallColijn(treeSym8, treeCat8), tolerance=1e-06)
  expect_equal(8L, KendallColijn(treeSym8, treeTac8), tolerance=1e-06)
  expect_equal(8L, KendallColijn(treeCat8, treeTac8), tolerance=1e-06)

  expect_equal(5.291503, KendallColijn(treeSym8, treeAb.Cdefgh), tolerance=1e-06)
  expect_equal(4.358899, KendallColijn(treeSym8, treeAbc.Defgh), tolerance=1e-06)
  expect_equal(5L, KendallColijn(treeSym8, treeAcd.Befgh), tolerance=1e-06)
  expect_equal(3.464102, KendallColijn(treeSym8, treeAbcd.Efgh), tolerance=1e-06)
  expect_equal(3L, KendallColijn(treeSym8, treeTwoSplits), tolerance=1e-06)
  expect_equal(2.828427, KendallColijn(treeAbc.Defgh, treeTwoSplits), tolerance=1e-06)
})
