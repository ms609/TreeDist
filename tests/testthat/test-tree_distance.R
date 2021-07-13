context('tree_distance.R')

# Labels in different order to confound as.Splits
treeSym8 <- ape::read.tree(text='((e, (f, (g, h))), (((a, b), c), d));')
treeBal8 <- ape::read.tree(text='(((e, f), (g, h)), ((a, b), (c, d)));')
treeOpp8 <- ape::read.tree(text='(((a, f), (c, h)), ((g, b), (e, d)));')
treesSBO8 <- structure(list(treeSym8, treeBal8, treeOpp8), 
                            class = 'multiPhylo')
treesSSBB8 <- structure(list(treeSym8, treeSym8, treeBal8, treeBal8), 
                            class = 'multiPhylo')

treeCat8 <- ape::read.tree(text='((((h, g), f), e), (d, (c, (b, a))));')
treeTac8 <- ape::read.tree(text='((((e, c), g), a), (h, (b, (d, f))));')
treeStar8 <- ape::read.tree(text='(e, c, g, h, b, a, d, f);')

treeAb.Cdefgh <- ape::read.tree(text='((a, b), (c, d, e, f, g, h));')
treeAbc.Defgh <- ape::read.tree(text='((a, b, c), (d, e, f, g, h));')
treeAcd.Befgh <- ape::read.tree(text='((a, c, d), (b, e, f, g, h));')
treeAbcd.Efgh <- ape::read.tree(text='((a, b, c, d), (e, f, g, h));')
treeTwoSplits <- ape::read.tree(text="(((a, b), c, d), (e, f, g, h));")

testTrees <- c(treesSBO8, treeCat8, treeTac8, treeStar8, treeAb.Cdefgh,
               treeAbc.Defgh, treeAbcd.Efgh, treeAcd.Befgh, treeTwoSplits)

test_that("Split compatibility is correctly established", {
  expect_true(SplitsCompatible(as.logical(c(0,0,1,1,0)), 
                               as.logical(c(0,0,1,1,0))))
  expect_true(SplitsCompatible( as.logical(c(0,0,1,1,0)), 
                               !as.logical(c(0,0,1,1,0))))
  expect_true(SplitsCompatible(as.logical(c(0,0,1,1,0)), 
                               as.logical(c(1,0,1,1,0))))
  expect_true(SplitsCompatible(!as.logical(c(0,0,1,1,0)), 
                                as.logical(c(1,0,1,1,0))))
  expect_false(SplitsCompatible(as.logical(c(0,0,1,1,0)), 
                                as.logical(c(1,1,0,1,0))))
})

methodsToTest <- list(
  SharedPhylogeneticInfo,
  DifferentPhylogeneticInfo,
  MatchingSplitInfo,
  MatchingSplitInfoDistance,
  MutualClusteringInfo,
  ClusteringInfoDistance,
  NyeSimilarity,
  JaccardRobinsonFoulds,
  MatchingSplitDistance,
  RobinsonFoulds,
  InfoRobinsonFoulds,
  KendallColijn # List last: requires rooted trees.
)

NormalizationTest <- function (FUNC, ...) {
  expect_equal(c(1L, 1L), 
               FUNC(treesSSBB8, normalize = TRUE, ...)[c(1, 6)],
               tolerance = 1e-7)
}

test_that('Bad labels cause error', {
  treeBadLabel8 <- ape::read.tree(text='((a, b, c, D), (e, f, g, h));')
  lapply(methodsToTest, function(Func) 
    expect_error(Func(treeSym8, treeBadLabel8)))
})

test_that('Size mismatch causes error', {
  treeSym7 <- ape::read.tree(text='((e, (f, g)), (((a, b), c), d));')
  splits7 <- as.Splits(treeSym7)
  splits8 <- as.Splits(treeSym8)
  
  lapply(methodsToTest, function(Func) 
    expect_error(Func(treeSym8, treeSym7)))
  
  lapply(methodsToTest, function(Func) 
    expect_error(Func(treeSym7, treeSym8)))
  
  expect_error(MeilaVariationOfInformation(splits7, splits8))
  
  Test <- function (Func) {
    expect_error(Func(splits8, as.Splits(BalancedTree(9)), 8))
  }
  Test(cpp_robinson_foulds_distance)
  Test(cpp_robinson_foulds_info)
  Test(cpp_matching_split_distance)
  Test(cpp_jaccard_similarity)
  Test(cpp_msi_distance)
  Test(cpp_mutual_clustering)
  Test(cpp_shared_phylo)
})

test_that('Metrics handle polytomies', {
  polytomy8 <- ape::read.tree(text='(a, b, c, d, e, f, g, h);')
  lapply(list(SharedPhylogeneticInfo, MutualClusteringInfo,
              MatchingSplitDistance, NyeSimilarity),
         function (Func) expect_equal(0, Func(treeSym8, polytomy8)))
})

#Func <- ClusteringInfoDistance # FUNC =
test_that('Output dimensions are correct', {
  list1 <- list(sym = treeSym8, bal = treeBal8)
  list2 <- list(sym = treeSym8, abc = treeAbc.Defgh, abcd = treeAbcd.Efgh)
  dimNames <- list(c('sym', 'bal'), c('sym', 'abc', 'abcd'))
  
  Test <- function (Func) {
    allPhylo <- matrix(
      c(Func(treeSym8, treeSym8),      Func(treeBal8, treeSym8),
        Func(treeSym8, treeAbc.Defgh), Func(treeBal8, treeAbc.Defgh),
        Func(treeSym8, treeAbcd.Efgh), Func(treeBal8, treeAbcd.Efgh)),
      2L, 3L, dimnames = dimNames)
    phylo1 <- matrix(c(Func(treeSym8, list2), Func(treeBal8, list2)),
                     byrow = TRUE, 2L, 3L, dimnames = dimNames)
    phylo2 <- matrix(c(Func(list1, treeSym8), Func(list1, treeAbc.Defgh),
                       Func(list1, treeAbcd.Efgh)), 2L, 3L, dimnames = dimNames)
    noPhylo <- Func(list1, list2)
    expect_equal(allPhylo, phylo1)
    expect_equal(allPhylo, phylo2)
    expect_equal(allPhylo, noPhylo)
  }
  
  lapply(methodsToTest, Test)
})

test_that('Robinson Foulds Distance is correctly calculated', {
  RFTest <- function (t1, t2) {
    expect_equal(suppressMessages(phangorn::RF.dist(t1, t2)),
                 RobinsonFoulds(t1, t2))
    
    expected <- RobinsonFoulds(t1, t2, reportMatching = TRUE, similarity = TRUE)
    attr(expected, 'pairScores') <- attr(expected, 'pairScores') == 0L
    expect_equal(expected, RobinsonFouldsMatching(t1, t2))
  }
  RFTest(treeSym8, treeSym8)
  RFTest(treeSym8, treeStar8)
  RFTest(treeStar8, treeStar8)
  RFTest(treeAb.Cdefgh, treeAbc.Defgh)
  RFTest(treeAb.Cdefgh, treeAbcd.Efgh)
  
  # at 2020-10, RF uses Day algorithm if tree2 = null; old algo if tree2 = tree1.
  expect_equivalent(RobinsonFoulds(testTrees, testTrees),
                    as.matrix(RobinsonFoulds(testTrees)))
  
  # Invariant to tree description order
  sq_pectinate <- ape::read.tree(text='((((((1, 2), 3), 4), 5), 6), (7, (8, (9, (10, 11)))));')
  shuffle1 <- ape::read.tree(text='(((((1, 5), 2), 6), (3, 4)), ((8, (7, 9)), (10, 11)));')
  shuffle2 <- ape::read.tree(text='(((8, (7, 9)), (10, 11)), ((((1, 5), 2), 6), (3, 4)));')
  RFTest(shuffle1, sq_pectinate)
  RFTest(sq_pectinate, shuffle1)
  RFTest(shuffle1, shuffle2)
  RFTest(shuffle1, sq_pectinate)
  RFTest(shuffle2, sq_pectinate)
})

test_that('Shared Phylogenetic Info is correctly calculated', {
  
  expect_equal(5.529821, tolerance = 1e-7,
               cpp_shared_phylo(
                 as.Splits(as.logical(c(1, 1, 1, 1, 0, 0, 0, 0))),
                 as.Splits(as.logical(c(1, 1, 1, 1, 0, 0, 0, 0))),
                 8L)$score)
  
  expect_equal(0.2895066, tolerance = 1e-7,
               cpp_shared_phylo(
                 as.Splits(as.logical(c(1, 1, 0, 0, 0, 0, 0, 0))),
                 as.Splits(as.logical(c(0, 0, 1, 1, 0, 0, 0, 0))),
                 8L)$score)
  
  expect_equal(1.137504, tolerance = 1e-6,
               cpp_shared_phylo(
                 as.Splits(as.logical(c(1, 1, 0, 0, 0, 0, 0, 0))),
                 as.Splits(as.logical(c(1, 1, 1, 1, 0, 0, 0, 0))),
                 8L)$score)
  
  expect_equal(3.45943, tolerance = 1e-6,
               cpp_shared_phylo(
                 as.Splits(as.logical(c(1, 1, 0, 0, 0, 0, 0, 0))),
                 as.Splits(as.logical(c(1, 1, 0, 0, 0, 0, 0, 0))),
                 8L)$score)
  
  expect_equal(22.53747, tolerance = 1e-05,
               SharedPhylogeneticInfo(treeSym8, treeSym8, normalize = FALSE))
  
  expect_equal(1, tolerance = 1e-05,
               SharedPhylogeneticInfo(treeSym8, treeSym8, normalize = TRUE))
  
  expect_equal(0,
               SharedPhylogeneticInfo(treeSym8, treeStar8, normalize = TRUE))

  expect_equal(0,
               SharedPhylogeneticInfo(treeStar8, treeStar8, normalize = FALSE))

  expect_equal(NaN, # Division by zero
               SharedPhylogeneticInfo(treeStar8, treeStar8, normalize = TRUE))
  
  expect_equal(13.75284, SharedPhylogeneticInfo(treeSym8, treeBal8), tolerance=1e-05)
  
  expect_equal(DifferentPhylogeneticInfo(treeSym8, treeAcd.Befgh),
               DifferentPhylogeneticInfo(treeAcd.Befgh, treeSym8), tolerance=1e-05)
  
  expect_equal(0, DifferentPhylogeneticInfo(treeSym8, treeSym8, normalize = TRUE))
  
  infoSymBal <- SplitwiseInfo(treeSym8) + SplitwiseInfo(treeBal8)
  expect_equal(infoSymBal - 13.75284 - 13.75284, tolerance = 1e-05,
    DifferentPhylogeneticInfo(treeSym8, treeBal8, normalize = TRUE) * infoSymBal)
  
  expect_equal(22.53747 + SharedPhylogeneticInfo(treeAcd.Befgh, treeAcd.Befgh) - 
                 (2 * SharedPhylogeneticInfo(treeSym8, treeAcd.Befgh)), 
               DifferentPhylogeneticInfo(treeSym8, treeAcd.Befgh), 
               tolerance=1e-06)
  
  expect_equal(-log2(945/10395), 
               SharedPhylogeneticInfo(treeSym8, treeAb.Cdefgh),
               tolerance = 1e-06)
  
  expect_equal(22.53747 + SharedPhylogeneticInfo(treeBal8, treeBal8) - 13.75284 - 13.75284, 
               DifferentPhylogeneticInfo(treeSym8, treeBal8), tolerance=1e-05)
  
  expect_equal(-log2(945/10395),
               SharedPhylogeneticInfo(treeSym8, treeAb.Cdefgh),
               tolerance = 1e-06)
  
  expect_equal(-log2(315/10395),
               SharedPhylogeneticInfo(treeSym8, treeAbc.Defgh),
               tolerance = 1e-06)
  
  expect_equal(0, DifferentPhylogeneticInfo(treeSym8, treeSym8))
  
  expect_equal(SplitwiseInfo(treeSym8) - SplitwiseInfo(treeAcd.Befgh),
               DifferentPhylogeneticInfo(treeSym8, treeAbc.Defgh),
               tolerance = 1e-06)
  
  
  # Test symmetry of small vs large splits
  expect_equal(SharedPhylogeneticInfo(treeSym8, treeAbc.Defgh),
               SharedPhylogeneticInfo(treeAbc.Defgh, treeSym8))
  expect_equal(-log2(225/10395), SharedPhylogeneticInfo(treeSym8, treeAbcd.Efgh))
  expect_equal(-log2(225/10395) - log2(945/10395),
               SharedPhylogeneticInfo(treeSym8, treeTwoSplits),
               tolerance = 1e-7)
  expect_equal(SplitSharedInformation(8, 4, 3),
               SharedPhylogeneticInfo(treeTwoSplits, treeAbc.Defgh),
               tolerance = 1e-7)
  expect_equal(SplitInformation(4, 4) + SplitInformation (3, 5) - 
               (2 * SplitSharedInformation(8, 4, 3)),
               SplitDifferentInformation(8, 4, 3),
               tolerance=1e-07)
  
  expect_equal(SharedPhylogeneticInfo(treeSym8, list(treeSym8, treeBal8)), 
               SharedPhylogeneticInfo(list(treeSym8, treeBal8), treeSym8),
               tolerance = 1e-7)
  
  # Test tree too large to cache
  set.seed(101)
  t1 <- ape::rtree(101)
  t2 <- ape::rtree(101, rooted = FALSE)
  expect_equal(SharedPhylogeneticInfo(t1, t2),
               SharedPhylogeneticInfo(t2, t1))
})

test_that('MatchingSplitInfo() is correctly calculated', {
  BinaryToSplit <- function (binary) matrix(as.logical(binary))
  expect_equal(log2(3),
               MatchingSplitInfoSplits(
                 as.Splits(c(rep(TRUE, 2), rep(FALSE, 6))),
                 as.Splits(c(FALSE, FALSE, rep(TRUE, 2), rep(FALSE, 4)))),
               tolerance = 1e-7)
  expect_equal(log2(3),
               MatchingSplitInfoSplits(
                 as.Splits(c(rep(FALSE, 6), rep(TRUE, 2))),
                 as.Splits(c(FALSE, FALSE, rep(TRUE, 2), rep(FALSE, 4)))),
               tolerance = 1e-7)
  expect_equal(log2(3), cpp_msi_distance(
    as.Splits(c(rep(TRUE, 2), rep(FALSE, 6))),
    as.Splits(c(FALSE, FALSE, rep(TRUE, 2), rep(FALSE, 4))),
    8L)$score, tolerance = 1e-7)
  expect_equal(log2(3), cpp_msi_distance(
    as.Splits(rep(c(FALSE, TRUE), each = 4L)),
    as.Splits(rep(c(FALSE, TRUE), 4L)),
    8L)$score, tolerance = 1e-7)
  
  
  expect_equal(SharedPhylogeneticInfo(treeSym8, treeSym8),
               MatchingSplitInfo(treeSym8, treeSym8), tolerance = 1e-05)
  
  expect_equal(0, MatchingSplitInfo(treeSym8, treeStar8))
  expect_equal(0, MatchingSplitInfo(treeStar8, treeStar8))
  
  expect_equal(MatchingSplitInfo(treeAb.Cdefgh, treeAbc.Defgh),
               MatchingSplitInfo(treeAbc.Defgh, treeAb.Cdefgh))
  expect_equal(MatchingSplitInfo(treeAbcd.Efgh, treeAb.Cdefgh),
               MatchingSplitInfo(treeAb.Cdefgh, treeAbcd.Efgh))
  expect_equal(-(TreeTools::Log2TreesMatchingSplit(2, 5) - Log2Unrooted.int(7)), 
               MatchingSplitInfo(treeAb.Cdefgh, treeAbc.Defgh),
               tolerance = 1e-7)
  expect_true(MatchingSplitInfo(treeSym8, treeBal8) > 
                MatchingSplitInfo(treeSym8, treeOpp8))
  expect_equal(0, MatchingSplitInfoDistance(treeSym8, treeSym8))
  NormalizationTest(MatchingSplitInfo)
})

test_that("Shared Phylogenetic Information is correctly estimated", {
  exp <- ExpectedVariation(treeSym8, treeAbc.Defgh, samples = 1000L)
  tol <- exp[, 'Std. Err.'] * 2
  # Expected values calculated with 100k samples
  expect_equal(1.175422, exp['SharedPhylogeneticInfo', 'Estimate'], 
               tolerance = tol[1])
  expect_equal(3.099776, exp['MatchingSplitInfo', 'Estimate'], 
               tolerance = tol[2])
  expect_equal(25.231023, exp['DifferentPhylogeneticInfo', 'Estimate'], 
               tolerance = tol[3])
  expect_equal(21.382314, exp['MatchingSplitInfoDistance', 'Estimate'], 
               tolerance = tol[4])
  expect_equal(exp[, 'sd'], exp[, 'Std. Err.'] * sqrt(exp[, 'n']))
})

test_that('Clustering information is correctly calculated', {
  expect_equal(Entropy(c(3, 5) / 8) * 2 - Entropy(c(0, 0, 3, 5) / 8),
               cpp_mutual_clustering(
                 as.Splits(as.logical(c(1, 1, 1, 0, 0, 0, 0, 0))),
                 as.Splits(as.logical(c(1, 1, 1, 0, 0, 0, 0, 0))),
                 8L)$score,
               tolerance = 1e-7)
  
  expect_equal(Entropy(c(2, 6) / 8) * 2 - Entropy(c(0, 2, 2, 4) / 8),
               cpp_mutual_clustering(
                 as.Splits(as.logical(c(1, 1, 0, 0, 0, 0, 0, 0))),
                 as.Splits(as.logical(c(0, 0, 1, 1, 0, 0, 0, 0))),
                 8L)$score, tolerance = 1e-7)
  
  expect_equal(Entropy(c(5, 4) / 9) + Entropy(c(3, 6) / 9) -
                      Entropy(c(3, 2, 0, 4) / 9),
               cpp_mutual_clustering(
                 as.Splits(as.logical(c(1, 1, 1, 1, 1, 0, 0, 0, 0))),
                 as.Splits(as.logical(c(0, 0, 1, 1, 1, 0, 0, 0, 0))),
                 9L)$score,
               tolerance = 1e-7)
  
  expect_equal(Entropy(c(4, 4) / 8) * 2 - Entropy(c(2, 2, 2, 2) / 8),
               cpp_mutual_clustering(
                 as.Splits(as.logical(c(1, 1, 1, 1, 0, 0, 0, 0))),
                 as.Splits(as.logical(c(1, 0, 1, 0, 1, 0, 1, 0))),
                 8L)$score,
               tolerance = 1e-7)
  
  expect_equal(Entropy(c(4, 4) / 8) * 2 - Entropy(c(0, 0, 4, 4) / 8),
               cpp_mutual_clustering(
                 as.Splits(as.logical(c(1, 1, 1, 1, 0, 0, 0, 0))),
                 as.Splits(as.logical(c(1, 1, 1, 1, 0, 0, 0, 0))),
                 8L)$score, 
               tolerance = 1e-7)
  
  expect_equal(ClusteringEntropy(treeSym8),
               MutualClusteringInfo(treeSym8, treeSym8),
               tolerance = 1e-05)
  expect_equal(8 * ClusteringEntropy(treeSym8), ClusteringInfo(treeSym8))
  expect_equal(0, MutualClusteringInfo(treeSym8, treeStar8))
  expect_equal(0, MutualClusteringInfo(treeStar8, treeStar8))
  
  expect_equal(TreeDistance(treeSym8, treeBal8),
               ClusteringInfoDistance(treeSym8, treeBal8, normalize = TRUE))
  expect_equal(1, MutualClusteringInfo(treeSym8, treeSym8, normalize = TRUE),
               tolerance = 1e-7)
  expect_true(MutualClusteringInfo(treeSym8, treeBal8, normalize = pmin) >
                MutualClusteringInfo(treeSym8, treeBal8, normalize = pmax))
  expect_equal(ClusteringEntropy(treeSym8) + ClusteringEntropy(treeBal8) -
                 (2 * MutualClusteringInfo(treeBal8, treeSym8)),
               ClusteringInfoDistance(treeSym8, treeBal8), tolerance = 1e-05)
  expect_equal(MutualClusteringInfo(treeAb.Cdefgh, treeAbc.Defgh),
               MutualClusteringInfo(treeAbc.Defgh, treeAb.Cdefgh),
               tolerance = 1e-05)
  
  
  
  # Different resolution
  randomBif20 <- structure(list(
    edge = structure(c(21L, 21L, 22L, 23L, 24L, 25L, 26L, 27L, 28L, 29L, 30L,
                       31L, 32L, 32L, 31L, 30L, 29L, 33L, 34L, 34L, 33L, 28L, 
                       35L, 36L, 36L, 35L, 27L, 26L, 37L, 37L, 25L, 38L, 38L,
                       39L, 39L, 24L, 23L, 22L, 1L, 22L, 23L, 24L, 25L, 26L, 
                       27L, 28L, 29L, 30L, 31L, 32L, 2L, 14L, 7L, 10L, 33L, 34L,
                       4L, 6L, 8L, 35L, 36L, 13L, 16L, 18L, 17L, 37L, 5L, 15L,
                       38L, 11L, 39L, 12L, 19L, 9L, 3L, 20L),
                     .Dim = c(38L, 2L)), Nnode = 19L,
    tip.label = c("t1", "t2", "t3", "t4", "t5", "t6", "t7", "t8", "t9", "t10",
                  "t11", "t12", "t13", "t14", "t15", "t16", "t17", "t18", "t19",
                  "t20"), br = NULL), class = "phylo")
  threeAwayPoly <- structure(
    list(edge = structure(c(21L, 22L, 23L, 24L, 25L, 26L, 27L, 28L, 29L, 29L, 
                            28L, 27L, 26L, 30L, 30L, 30L, 26L, 31L, 31L, 25L, 
                            32L, 33L, 33L, 32L, 25L, 25L, 24L, 34L, 34L, 34L,
                            23L, 22L, 21L, 22L, 23L, 24L, 25L, 26L, 27L, 28L,
                            29L, 2L, 8L, 14L, 10L, 30L, 13L, 16L, 18L, 31L, 4L,
                            6L, 32L, 33L, 15L, 20L, 5L, 7L, 17L, 34L, 11L, 12L,
                            19L, 9L, 3L, 1L), .Dim = c(33L, 2L)),
         tip.label = c("t1", "t2", "t3", "t4", "t5", "t6", "t7", "t8", "t9",
                       "t10", "t11", "t12", "t13", "t14", "t15", "t16", "t17",
                       "t18", "t19", "t20"),
         Nnode = 14L), class = "phylo")
  expect_equal(
    MutualClusteringInfo(threeAwayPoly, randomBif20),
    MutualClusteringInfo(randomBif20, threeAwayPoly))
  match <- MutualClusteringInfo(randomBif20, threeAwayPoly, reportMatching = TRUE)
  expect_equal(c(NA, NA,  1,  2, NA,  3,  7, 11, 10,  4,  6,  9,  8, NA,  5, 12, NA),
               attr(match, 'matching'))
  
  library('TreeTools')
  expect_equal(ClusteringEntropy(BalancedTree(64)),
               MutualClusteringInfo(BalancedTree(64), BalancedTree(64)))
  expect_equal(ClusteringEntropy(BalancedTree(644)),
               MutualClusteringInfo(BalancedTree(644), BalancedTree(644)))
  
  expect_gt(ClusteringEntropy(BalancedTree(64)),
            MutualClusteringInfo(BalancedTree(64), PectinateTree(64)))
  expect_gt(ClusteringEntropy(BalancedTree(644)),
            MutualClusteringInfo(BalancedTree(644), PectinateTree(644)))
  
  NormalizationTest(MutualClusteringInfo)
})

test_that("Matchings are correct", {
  
  # Different resolution: used to cause memory leak
  randomBif20 <- structure(list(
    edge = structure(c(21L, 21L, 22L, 23L, 24L, 25L, 26L, 27L, 28L, 29L, 30L,
                       31L, 32L, 32L, 31L, 30L, 29L, 33L, 34L, 34L, 33L, 28L, 
                       35L, 36L, 36L, 35L, 27L, 26L, 37L, 37L, 25L, 38L, 38L,
                       39L, 39L, 24L, 23L, 22L, 1L, 22L, 23L, 24L, 25L, 26L, 
                       27L, 28L, 29L, 30L, 31L, 32L, 2L, 14L, 7L, 10L, 33L, 34L,
                       4L, 6L, 8L, 35L, 36L, 13L, 16L, 18L, 17L, 37L, 5L, 15L,
                       38L, 11L, 39L, 12L, 19L, 9L, 3L, 20L),
                     .Dim = c(38L, 2L)), Nnode = 19L,
    tip.label = c("t1", "t2", "t3", "t4", "t5", "t6", "t7", "t8", "t9", "t10",
                  "t11", "t12", "t13", "t14", "t15", "t16", "t17", "t18", "t19",
                  "t20"), br = NULL), class = "phylo")
  threeAwayPoly <- structure(
    list(edge = structure(c(21L, 22L, 23L, 24L, 25L, 26L, 27L, 28L, 29L, 29L, 
                            28L, 27L, 26L, 30L, 30L, 30L, 26L, 31L, 31L, 25L, 
                            32L, 33L, 33L, 32L, 25L, 25L, 24L, 34L, 34L, 34L,
                            23L, 22L, 21L, 22L, 23L, 24L, 25L, 26L, 27L, 28L,
                            29L, 2L, 8L, 14L, 10L, 30L, 13L, 16L, 18L, 31L, 4L,
                            6L, 32L, 33L, 15L, 20L, 5L, 7L, 17L, 34L, 11L, 12L,
                            19L, 9L, 3L, 1L), .Dim = c(33L, 2L)),
         tip.label = c("t1", "t2", "t3", "t4", "t5", "t6", "t7", "t8", "t9",
                       "t10", "t11", "t12", "t13", "t14", "t15", "t16", "t17",
                       "t18", "t19", "t20"),
         Nnode = 14L), class = "phylo")
  
  expect_equal(
    MutualClusteringInfo(threeAwayPoly, randomBif20),
    MutualClusteringInfo(randomBif20, threeAwayPoly))
  
  
  t1 <- PectinateTree(letters[1:11])
  t2 <- ape::read.tree(text = '(a, (c, (b, (d, e, ((g, h, f), (k, (j, i)))))));')
  t3 <- CollapseNode(PectinateTree(c(letters[11], letters[1:10])), 16:19)
  s1 <- as.Splits(t1)
  s2 <- as.Splits(t2, t1)
  s3 <- as.Splits(t3, t1)
  n1 <- dim(s1)[1]
  n2 <- dim(s2)[1]
  n3 <- dim(s3)[1]
  n <- NTip(s1)
  
  # Plot
  # par(mfrow = 2:1, cex = 0.9, mar = rep(0,4))
  # JRF2T <- function(...) JaccardRobinsonFoulds(..., k = 2)
  # JRF2F <- function(...) JaccardRobinsonFoulds(..., k = 2, allowConflict = FALSE)
  # VisualizeMatching(MatchingSplitDistance, t1, t2, setPar=F)
  # LabelSplits(t2, setNames(1:6, names(s2)), adj = 2)
  # VisualizeMatching(MatchingSplitDistance, t2, t1, setPar=F)
  # LabelSplits(t1, setNames(1:8, names(s1)), adj = 2)
  
  
  Test <- function (CppFn, x12, x21, ...) {
    
    r12 <- CppFn(s1, s2, n, ...)
    r21 <- CppFn(s2, s1, n, ...)
    r13 <- CppFn(s1, s3, n, ...)
    r31 <- CppFn(s3, s1, n, ...)
    expect_equal(r12$score, r21$score)
    expect_equal(r13$score, r31$score)
    
    m12 <- r12$matching
    m21 <- r21$matching
    
    expect_equal(n1, length(m12))
    expect_equal(length(m12[!is.na(m12)]), length(unique(m12[!is.na(m12)])))
    expect_equal(n2, length(m21))
    expect_equal(length(m21[!is.na(m21)]), length(unique(m21[!is.na(m21)])))
    expect_lte(dim(s1)[1] - dim(s2)[1], sum(is.na(m12)))
    
    
    m13 <- r13$matching
    m31 <- r31$matching
    expect_equal(n1, length(m13))
    expect_equal(length(m13[!is.na(m13)]), length(unique(m13[!is.na(m13)])))
    expect_equal(n3, length(m31))
    expect_equal(length(m31[!is.na(m31)]), length(unique(m31[!is.na(m31)])))
    expect_lte(dim(s1)[1] - dim(s3)[1], sum(is.na(m13)))
    
    for (i in seq_along(m12)) expect_true(m12[i] %in% x12[[i]])
    for (i in seq_along(m21)) expect_true(m21[i] %in% x21[[i]])
    
  }
  
  Test(TreeDist:::cpp_robinson_foulds_distance,
       list(NA, 2, NA, 3, NA, NA, 5, NA),
       list(NA, 2, 4, NA, 7, NA)
       )
  Test(TreeDist:::cpp_robinson_foulds_info,
       list(NA, 2, NA, 3, NA, NA, 5, NA),
       list(NA, 2, 4, NA, 7, NA)
       )
  Test(TreeDist:::cpp_matching_split_distance,
       list(1, 2, 4, 3, NA, NA, 5, 6),
       list(1, 2, 5, 4, 7, 6)
       )
  Test(TreeDist:::cpp_jaccard_similarity,
       list(NA, 2, 1, 3, 4, 6, 5, NA),
       list(3, 2, 4, 5, 7, 6),
       k = 2,
       allowConflict = TRUE)
  Test(TreeDist:::cpp_jaccard_similarity,
       list(NA, 2, 1, 3, NA, 6, 5, 4),
       list(3, 2, 4, 1, 7, 6),
       k = 2,
       allowConflict = FALSE)
  Test(TreeDist:::cpp_msi_distance,
       list(NA, 2, 1, 4, 3, 6, 5, NA),
       list(3, 2, c(4, 5), c(4, 5), c(6, 7), c(7, 6))
       )
  Test(TreeDist:::cpp_shared_phylo,
       list(NA, 2, 4, 3, 1, 6, 5, NA),
       list(5, 2, 4, 3, 7, 6)
       )
  Test(TreeDist:::cpp_mutual_clustering, 
       list(4, 2, NA, 3, 6, NA, 5, 1),
       list(8, 2, 4, 5, 7, 1)
       )

})

test_that('Matching Split Distance is correctly calculated', {
  expect_equal(0L, MatchingSplitDistance(treeSym8, treeSym8))
  expect_equal(0L, MatchingSplitDistance(treeStar8, treeSym8))
  expect_equal(0L, MatchingSplitDistance(treeStar8, treeStar8))
  match0 <- MatchingSplitDistance(treeStar8, treeStar8, reportMatching = TRUE)
  expect_equivalent(rep(0L, 4), 
                    c(match0, vapply(attributes(match0), length, 0)))
  expect_equal(1L, MatchingSplitDistance(treeAb.Cdefgh, treeAbc.Defgh))
  expect_equal(2L, MatchingSplitDistance(treeAb.Cdefgh, treeAbcd.Efgh))
  
  splitAB <- as.Splits(c(rep(TRUE, 2), rep(FALSE, 7)))
  splitABC <- as.Splits(c(rep(TRUE, 3), rep(FALSE, 6)))
  splitAEF <- as.Splits(c(TRUE, rep(FALSE, 3), TRUE, TRUE, rep(FALSE, 3)))
  splitABCD <- as.Splits(c(rep(TRUE, 4), rep(FALSE, 5)))
  splitABCDE <- as.Splits(c(rep(TRUE, 5), rep(FALSE, 4)))
  splitAI <- as.Splits(c(TRUE, rep(FALSE, 7), TRUE))
  
  expect_equal(2L, MatchingSplitDistanceSplits(splitAB, splitAI))
  expect_equal(2L, MatchingSplitDistanceSplits(splitAB, splitABCD))
  expect_equal(3L, MatchingSplitDistanceSplits(splitAB, splitABCDE))
  expect_equal(4L, MatchingSplitDistanceSplits(splitABC, splitAEF))
  expect_equal(MatchingSplitDistanceSplits(splitABC, splitAEF),
               MatchingSplitDistanceSplits(splitAEF, splitABC))
  
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

test_that('NyeSimilarity is correctly calculated, and matches JRF', {
  listBalSym <- list(treeBal8, treeSym8)
  
  JRF <- function (..., sim = TRUE)
    JaccardRobinsonFoulds(..., k = 1, similarity = sim, allowConflict = TRUE)
  expect_equal(5L, NyeSimilarity(as.Splits(treeSym8), treeSym8))
  
  expect_equal(1, NyeSimilarity(treeSym8, treeSym8, normalize = TRUE))
  expect_equal(1, JRF(treeSym8, treeSym8, normalize = TRUE))
  
  expect_equal(0, NyeSimilarity(treeSym8, treeStar8, normalize = FALSE))
  expect_equal(0, NyeSimilarity(treeSym8, treeStar8, normalize = TRUE))
  expect_equal(0, JRF(treeSym8, treeStar8, normalize = TRUE))
  
  expect_equal(0, NyeSimilarity(treeStar8, treeStar8, normalize = FALSE))
  expect_equal(NaN, NyeSimilarity(treeStar8, treeStar8, normalize = TRUE, 
                                      normalizeMax = FALSE))
  
  expect_equal(c(3.8, 5), NyeSimilarity(treeSym8, listBalSym))
  expect_equal(2 / 3, NyeSimilarity(treeAb.Cdefgh, treeAbc.Defgh), 
               tolerance = 1e-7)
  expect_equal(2 * (1 / 3), tolerance = 1e-7, 
               NyeSimilarity(treeAb.Cdefgh, treeAbc.Defgh, similarity = FALSE))
  
  expect_equal(1L, NyeSimilarity(treeSym8, treeAbcd.Efgh, normalize = FALSE))
  expect_equal(1L / 5L, NyeSimilarity(treeSym8, treeAbcd.Efgh, normalize = 5L))
  expect_equal(0.2, JRF(treeSym8, treeAbcd.Efgh, normalize = 5L * 2L))
  expect_equal(1/3, NyeSimilarity(treeSym8, treeAbcd.Efgh, normalize = TRUE))
  expect_equal(1/3, JRF(treeSym8, treeAbcd.Efgh, normalize = TRUE))
  expect_equal(2/3, NyeSimilarity(treeSym8, treeAbcd.Efgh, similarity = FALSE,
                                  normalize = TRUE))
  expect_equal(2/3, JRF(treeSym8, treeAbcd.Efgh, sim = FALSE, normalize = TRUE))
  
  expect_equal(1L / ((5L + 1L) / 2L),
               NyeSimilarity(treeSym8, treeAbcd.Efgh, normalize = TRUE))
  expect_true(NyeSimilarity(treeSym8, treeBal8) > 
                NyeSimilarity(treeSym8, treeOpp8))
  
  NormalizationTest(NyeSimilarity)
})


test_that('Jaccard RF extremes tend to equivalent functions', {
  expect_equal(JaccardRobinsonFoulds(treeSym8, list(treeBal8, treeSym8),
                                     similarity = TRUE, k = 1L,
                                     allowConflict = TRUE),
               NyeSimilarity(treeSym8, list(treeBal8, treeSym8)) * 2L)
  
  expect_equal(JaccardRobinsonFoulds(treeSym8, list(treeBal8, treeSym8),
                                     similarity = FALSE, k = Inf),
               RobinsonFoulds(treeSym8, list(treeBal8, treeSym8)))
  
  expect_equal(JaccardRobinsonFoulds(treeSym8, list(treeBal8, treeSym8),
                                     similarity = FALSE, k = 999999),
               RobinsonFoulds(treeSym8, list(treeBal8, treeSym8)))
})

test_that('Jaccard RF is correctly calculated', {
  expect_equal(5L * 2L, JaccardRobinsonFoulds(treeSym8, treeSym8,
                                         k = 2, similarity = TRUE))
  expect_equal(c(3.32, 5) * 2L, 
               JaccardRobinsonFoulds(treeSym8, list(treeBal8, treeSym8),
                                     similarity = TRUE, k = 2))
  expect_equal(2 * 2, 3 * JaccardRobinsonFoulds(treeAb.Cdefgh, treeAbc.Defgh,
                                            similarity = TRUE),
               tolerance = 1e-7)
  expect_equal(1, JaccardRobinsonFoulds(treeSym8, treeSym8,
                                        similarity = TRUE, normalize = TRUE))
  expect_equal(0, JaccardRobinsonFoulds(treeSym8, treeSym8,
                                        similarity = FALSE, normalize = TRUE))
  expect_equal(1L * 2L, 
               JaccardRobinsonFoulds(treeSym8, treeAbcd.Efgh, similarity = TRUE,
                                     normalize = FALSE, k = 2))
  expect_equal(1L * 2L / 6L, 
               JaccardRobinsonFoulds(treeSym8, treeAbcd.Efgh, similarity = TRUE,
                                     normalize = TRUE, k = 2))
  expect_lt(JaccardRobinsonFoulds(treeSym8, treeBal8, k = 2),
            JaccardRobinsonFoulds(treeSym8, treeOpp8, k = 2))
  expect_lt(JaccardRobinsonFoulds(treeSym8, treeBal8, k = 3L),
            JaccardRobinsonFoulds(treeSym8, treeBal8, k = 4L))
  expect_lt(JaccardRobinsonFoulds(treeCat8, treeTac8, allowConflict = TRUE),
            JaccardRobinsonFoulds(treeCat8, treeTac8, allowConflict = FALSE))
  
  expect_equal(0, JaccardRobinsonFoulds(BalancedTree(64), BalancedTree(64)))
  expect_lt(0, JaccardRobinsonFoulds(BalancedTree(64), PectinateTree(64)))
  expect_equal(0, JaccardRobinsonFoulds(BalancedTree(264), BalancedTree(264)))
  expect_lt(0, JaccardRobinsonFoulds(BalancedTree(264), PectinateTree(264)))
})

test_that('RobinsonFoulds() is correctly calculated', {
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
  
  RFNtipTest <- function (nTip) {
    backLeaves <- paste0('t', rev(seq_len(nTip)))
    RFTest(TreeTools::PectinateTree(backLeaves), 
           TreeTools::BalancedTree(nTip))
  }
  RFNtipTest(10)
  RFNtipTest(32)
  RFNtipTest(50)
  RFNtipTest(64)
  RFNtipTest(67)
  RFNtipTest(128)
  RFNtipTest(1024)
  RFNtipTest(1027)
  
  NormalizationTest(RobinsonFoulds, similarity = TRUE)
  #TODO we may wish to revise this test once we implement diag = TRUE to 
  #allow similarities to be calculated on the diagonal.
  expect_equal(numeric(0), RobinsonFoulds(treeSym8, normalize = TRUE))
})


test_that('Robinson Foulds Info is correctly calculated', {
  expect_equal(22.53747 * 2L, tolerance = 1e-05,
               InfoRobinsonFoulds(treeSym8, treeSym8, similarity = TRUE,
                                  normalize = FALSE))
  expect_equal(0, tolerance = 1e-05,
               InfoRobinsonFoulds(treeSym8, treeSym8, normalize = TRUE))
  expect_equal(1, tolerance = 1e-05,
               InfoRobinsonFoulds(treeSym8, treeSym8, similarity = TRUE, 
                                  normalize = TRUE))
  expect_equal(24.9, tolerance = 0.01, 
               InfoRobinsonFoulds(treeSym8, treeBal8, similarity = TRUE))
  expect_equal(SplitwiseInfo(treeSym8) + SplitwiseInfo(treeBal8) -
                 InfoRobinsonFoulds(treeSym8, treeBal8, similarity = FALSE),
               InfoRobinsonFoulds(treeSym8, treeBal8, similarity = TRUE))
  expect_equal(-log2(945/10395) * 2,
               InfoRobinsonFoulds(treeSym8, treeAb.Cdefgh, similarity = TRUE))
  expect_equal(-log2(945/10395) * 2, 
               InfoRobinsonFoulds(treeSym8, treeAb.Cdefgh, similarity = TRUE))
  expect_equal(-log2(315/10395) * 2, 
               InfoRobinsonFoulds(treeSym8, treeAbc.Defgh, similarity = TRUE))
  
  # Test symmetry of small vs large splits
  expect_equal(InfoRobinsonFoulds(treeSym8, treeAbc.Defgh),
               InfoRobinsonFoulds(treeAbc.Defgh, treeSym8))
  expect_equal(-log2(225/10395) * 2, 
               InfoRobinsonFoulds(treeSym8, treeAbcd.Efgh, similarity = TRUE))
  expect_equal((-log2(225/10395) - log2(945/10395)) * 2,
               InfoRobinsonFoulds(treeSym8, treeTwoSplits, similarity = TRUE))
  expect_equal(InfoRobinsonFoulds(treeSym8, list(treeSym8, treeBal8)), 
               RobinsonFouldsInfo(list(treeSym8, treeBal8), treeSym8))
  
  # Check that large trees work
  expect_equal(0, InfoRobinsonFoulds(BalancedTree(64), BalancedTree(64)))
  expect_lt(0, InfoRobinsonFoulds(BalancedTree(64), PectinateTree(64)))
  expect_equal(0, InfoRobinsonFoulds(BalancedTree(129), BalancedTree(129)))
  expect_lt(0, InfoRobinsonFoulds(BalancedTree(129), PectinateTree(129)))
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

test_that('Multiple comparisons are correctly ordered', {
  nTrees <- 6L
  nTip <- 16L
  
  set.seed(0)
  trees <- lapply(rep(nTip, nTrees), ape::rtree, br=NULL)
  trees[[1]] <- TreeTools::BalancedTree(nTip)
  trees[[nTrees - 1L]] <- TreeTools::PectinateTree(nTip)
  class(trees) <- 'multiPhylo'
  
  expect_equivalent(phangorn::RF.dist(trees),
                    RobinsonFoulds(trees))
  
  # Test CompareAll
  expect_equivalent(as.matrix(phangorn::RF.dist(trees)),
                    as.matrix(CompareAll(trees, phangorn::RF.dist, 0L)))
  
  NNILoose <- function (x, y) NNIDist(x, y)['loose_upper']
  expect_equivalent(CompareAll(trees, NNILoose),
                    CompareAll(trees, NNIDist)$loose_upper)
})

test_that('Normalization occurs as documented', {
  library('TreeTools')
  tree1 <- BalancedTree(8)
  tree2 <- CollapseNode(PectinateTree(8), 12:13)
  
  info1 <- SplitwiseInfo(tree1) # 19.367
  info2 <- SplitwiseInfo(tree2) # 11.963
  
  ent1 <- ClusteringEntropy(tree1) # 4.245
  ent2 <- ClusteringEntropy(tree2) # 2.577
  
  # Phylogenetic information
  spi <- SharedPhylogeneticInfo(tree1, tree2, normalize = FALSE) # 9.64
  dpi <- DifferentPhylogeneticInfo(tree1, tree2, normalize = FALSE) # 12.04
  expect_equal(spi + spi + dpi, info1 + info2)
  
  expect_equal(SharedPhylogeneticInfo(tree1, tree2, normalize = TRUE),
               (spi + spi) / (info1 + info2))
  expect_equal(PhylogeneticInfoDistance(tree1, tree2, normalize = TRUE),
               dpi / (info1 + info2))
  
  # Matching split information
  mmsi <- MatchingSplitInfo(tree1, tree2, normalize = FALSE)
  msid <- MatchingSplitInfoDistance(tree1, tree2, normalize = FALSE)
  expect_equal(mmsi + mmsi + msid, info1 + info2)
  
  expect_equal(MatchingSplitInfo(tree1, tree2, normalize = TRUE),
               (mmsi + mmsi) / (info1 + info2))
  expect_equal(MatchingSplitInfoDistance(tree1, tree2, normalize = TRUE),
               msid / (info1 + info2))
  
  
  # Clustering information
  mci <- MutualClusteringInfo(tree1, tree2, normalize = FALSE)
  cid <- ClusteringInfoDistance(tree1, tree2, normalize = FALSE)
  expect_equal(mci + mci + cid, ent1 + ent2)
  
  expect_equal(MutualClusteringInfo(tree1, tree2, normalize = TRUE),
               (mci + mci) / (ent1 + ent2))
  expect_equal(ClusteringInfoDistance(tree1, tree2, normalize = TRUE),
               cid / (ent1 + ent2))
  
})

test_that("Independent of root position", {
  
  library('TreeTools')
  
  bal8 <- BalancedTree(8)
  pec8 <- PectinateTree(8)
  
  trees <- lapply(list(bal8, RootTree(bal8, 't4'),
                       pec8, RootTree(pec8, 't4')), UnrootTree)
  
  lapply(methodsToTest[-length(methodsToTest)], function (Method) {
    dists <- as.matrix(Method(trees))
    expect_equal(dists[1, 1], dists[1, 2])
    expect_equal(dists[1, 3], dists[1, 4])
    expect_equal(dists[1, 3], dists[2, 4])
    expect_equal(dists[2, 3], dists[2, 4])
    expect_equal(dists[3, 3], dists[3, 4])
  })

  
  Test <- function(Method, score = 0L, ...) {
    expect_equal(score, Method(trees[[1]], trees[[1]], ...))
    expect_equal(score, Method(trees[[1]], trees[[2]], ...))
    expect_equal(score, Method(trees[[3]], trees[[3]], ...))
  }
  
  Test(MASTSize, 8L, rooted = FALSE)
  # Tested further for NNIDist in test-tree_distance_nni.R
  Test(NNIDist, c(lower = 0, best_lower = 0, tight_upper = 0, best_upper = 0,
                  loose_upper = 0, fack_upper = 0, li_upper = 0))
  Test(SPRDist, c(spr = 0))
  
})
