library('TreeTools')
context('tree_distance_utilities.R')

test_that('Tree normalization works', {
  expect_equal(0.5, NormalizeInfo(5, 3, 5, how=TRUE, 
                                  InfoInTree=function(x, y) x + y, y = 1L))
  expect_equal(0:4 / c(1, 2, 3, 3, 3), 
               NormalizeInfo(0:4, 1:5, 3, InfoInTree = I, Combine = min))
  expect_equal(3/4, NormalizeInfo(unnormalized = 3, 1, 1, how = 4L))
  expect_equal(matrix(0:24 / 10, 5, 5), 
               NormalizeInfo(0:24, rep(5, 5), rep(5, 5), InfoInTree = I))
  expect_equal(matrix(1:20 / 10, 5, 4) * c(10,10,1,1,1, 10,10,1,1,1, rep(1,10)),
               NormalizeInfo(1:20, c(1, 1, 10, 10, 10), c(1, 1, 10, 10), 
                             InfoInTree = I, Combine = pmax))
})


test_that('CalculateTreeDistance handles splits appropriately', {
  set.seed(101)
  tree10 <- ape::rtree(10)
  tree10.1 <- ape::rtree(10)
  splits10 <- as.Splits(tree10)
  splits10.1 <- as.Splits(tree10.1)
  trees10.3 <- list(BalancedTree(tree10),
                    PectinateTree(tree10),
                    ape::rtree(10, br = NULL))
  trees10.4 <- list(PectinateTree(tree10),
                    BalancedTree(tree10),
                    ape::rtree(10, br = NULL),
                    ape::rtree(10, br = NULL))
  splits10.3 <- as.Splits(trees10.3)
  splits10.4 <- as.Splits(trees10.4)
  
  expect_equal(0, CalculateTreeDistance(RobinsonFouldsSplits, splits10, splits10))
  expect_equal(seq_len(7),
               attr(CalculateTreeDistance(RobinsonFouldsSplits, 
                                          splits10, splits10, TRUE), 'matching')
               )
  expect_equal(CalculateTreeDistance(RobinsonFouldsSplits, splits10, splits10.1),
               CalculateTreeDistance(RobinsonFouldsSplits, splits10.1, splits10))
  expect_equal(CalculateTreeDistance(RobinsonFouldsSplits, splits10, splits10.3),
               CalculateTreeDistance(RobinsonFouldsSplits, splits10.3, splits10))
  
  expect_equivalent(
    CalculateTreeDistance(RobinsonFouldsSplits, splits10.3, splits10.3),
    CalculateTreeDistance(RobinsonFouldsSplits, splits10.3, trees10.3))
  
  expect_equivalent(
    CalculateTreeDistance(RobinsonFouldsSplits, trees10.3, splits10.3),
    CalculateTreeDistance(RobinsonFouldsSplits, trees10.3, trees10.3))
  
  expect_equivalent(
    CalculateTreeDistance(RobinsonFouldsSplits, trees10.3, trees10.3)[, c(1, 3, 2)],
    CalculateTreeDistance(RobinsonFouldsSplits, trees10.3, trees10.3[c(1, 3, 2)]),
  )
  
  mat <- matrix(NA, 3, 4)
  for (i in 1:3) for (j in 1:4) {
    mat[i, j] <- CalculateTreeDistance(RobinsonFouldsSplits, splits10.3[[i]], splits10.4[[j]])
  }
  
  expect_equivalent(mat, CalculateTreeDistance(RobinsonFouldsSplits, splits10.3, splits10.4))
  
  expect_equivalent(
    CalculateTreeDistance(RobinsonFouldsSplits, splits10.4, splits10.3),
    t(CalculateTreeDistance(RobinsonFouldsSplits, splits10.3, splits10.4)))
})

test_that('Matches are reported', {
  # Trees copied from test-tree_distance.R
  treeSym8 <- ape::read.tree(text='((e, (f, (g, h))), (((a, b), c), d));')
  treeBal8 <- ape::read.tree(text='(((e, f), (g, h)), ((a, b), (c, d)));')
  treeTwoSplits <- ape::read.tree(text="(((a, b), c, d), (e, f, g, h));")
  
  A <- TRUE
  B <- FALSE
  splitsA <- as.Splits(c(A, A, B, B, B, B, B, B), tipLabels = letters[1:8])
  splitsB <- as.Splits(matrix(c(A, A, A, A, A, A, B, B, 
                          A, A, A, B, B, B, B, B,
                          B, B, A, A, A, A, A, A,
                          A, A, A, A, A, B, A, B), nrow=4, byrow=TRUE),
                       tipLabels = letters[1:8])
  
  expect_equal("a b | c d e f g h => a b c | d e f g h",
               ReportMatching(splitsA, splitsB[[2]]))
  
  splits1 <- as.Splits(treeSym8)
  splits2 <- as.Splits(treeBal8)
  
  matchedSplits <- match.Splits(splits1, splits2)
  cs <- CompatibleSplits(splits1, splits2)
  cs[, matchedSplits] <- FALSE
  unmatched <- is.na(matchedSplits)
  matchedSplits[unmatched] <- apply(cs[unmatched, ], 1, which)
  expect_equal(matchedSplits, attr(GeneralizedRF(as.Splits(treeSym8),
                                                    as.Splits(treeBal8), 8L, 
                                                    cpp_shared_phylo,
                                                    maximize = TRUE, 
                                                    reportMatching = TRUE),
                                      'matching'))
    
  expect_equal(matchedSplits, attr(
    SharedPhylogeneticInfoSplits(as.Splits(treeSym8), as.Splits(treeBal8),
                               reportMatching = TRUE),
    'matching'))
  
  Test <- function (Func, relaxed = FALSE, ...) {
    tree1 <- PectinateTree(letters[1:8])
    tree2 <- BalancedTree(letters[8:1])
    at <- attributes(Func(tree1, tree2, reportMatching = TRUE, ...))
    expect_equal(3L, length(at))
    
    splits1 <- as.Splits(tree1)
    splits2 <- as.Splits(tree2, tree1)
    
    matchedSplits <- match.Splits(splits1, splits2)
    if (relaxed) {
      expect_equal(matchedSplits[!is.na(matchedSplits)],
                   as.integer(at$matching[c(1, 3, 5)]))
    } else {
      cs <- CompatibleSplits(splits1, splits2)
      cs[, matchedSplits] <- FALSE
      unmatched <- is.na(matchedSplits)
      matchedSplits[unmatched] <- apply(cs[unmatched, ], 1, which)
      
      expect_equal(matchedSplits, as.integer(at$matching))
    }
    ghSplit <- at$matchedSplits[
      match.Splits(as.Splits(c(rep(FALSE, 6), TRUE, TRUE), letters[1:8]),
                   splits1[[which(!is.na(matchedSplits))]])]
    expect_equal('g h | a b c d e f => g h | a b c d e f', ghSplit)
    
    at <- attributes(Func(treeSym8, treeTwoSplits, reportMatching = TRUE, ...))
    expect_equal(3L, length(at))
    expect_equal(match.Splits(as.Splits(treeSym8), as.Splits(treeTwoSplits, treeSym8)),
                 as.integer(at$matching))
    expect_equal('a b | e f g h c d => a b | e f g h c d', at$matchedSplits[2])
  }
  
  Test(SharedPhylogeneticInfo)
  Test(DifferentPhylogeneticInfo)
  Test(MatchingSplitInfo)
  Test(MatchingSplitInfoDistance)
  Test(MutualClusteringInfo)
  Test(ClusteringInfoDistance)
  
  Test(NyeTreeSimilarity)
  Test(MatchingSplitDistance)
  Test(JaccardRobinsonFoulds, k = 2, arboreal = FALSE)
  Test(JaccardRobinsonFoulds, k = 2, arboreal = TRUE)
  Test(RobinsonFoulds, relaxed = TRUE)
  Test(InfoRobinsonFoulds, relaxed = TRUE)

    # Matching Split Distance matches differently:  
  at <- attributes(MatchingSplitDistance(treeSym8, treeBal8, 
                                         reportMatching = TRUE))
  expect_equal(3L, length(at))
  expect_equal(c(1:3, 5:4), as.integer(at$matching))
  expect_equal('a b | e f g h c d => a b | e f g h c d', at$matchedSplits[5])
  
  # Zero match:
  expect_equal('a b | c d .. a c | b d', 
               attr(SharedPhylogeneticInfo( 
                      ape::read.tree(text="((a, b), (c, d));"),
                      ape::read.tree(text="((a, c), (b, d));"), 
                      reportMatching = TRUE),
                    'matchedSplits'))
})

test_that('Matchings are calculated in both directions', {
  tree1 <- ape::read.tree(text='((1, 2), ((3, 4, (5, 9)), (6, (7, 8))));')
  tree2 <- ape::read.tree(text='((1, 2), ((3, (4, 5)), (6, (7, (8, 9)))));')
  splits1 <- as.Splits(tree1)
  splits2 <- as.Splits(tree2)
  nMatches <- min(length(splits1), length(splits2))
  
  Test <- function (Func, ...) {
    matching12 <- Func(tree1, tree2, reportMatching = TRUE, ...)
    expect_equal(nMatches, length(attr(matching12, 'matchedSplits')))
    
    matching21 <- Func(tree2, tree1, reportMatching = TRUE, ...)
    expect_equal(nMatches, length(attr(matching21, 'matchedSplits')))
  }
  
  Test(SharedPhylogeneticInfo)
  Test(DifferentPhylogeneticInfo)
  Test(MatchingSplitInfo)
  Test(MatchingSplitInfoDistance)
  Test(MutualClusteringInfo)
  Test(ClusteringInfoDistance)
  
  Test(NyeTreeSimilarity)
  Test(JaccardRobinsonFoulds, k = 2, arboreal = FALSE)
  Test(JaccardRobinsonFoulds, k = 2, arboreal = TRUE)
  Test(MatchingSplitDistance)
  
  nMatches <- 1L
  Test(RobinsonFoulds)
  Test(InfoRobinsonFoulds)
  
})

test_that('.TreeDistance supports all sizes', {
  expect_equal(rbind(
    bal = MASTSize(BalancedTree(7), as.phylo(0:3, 7)),
    pec = MASTSize(as.phylo(0:3, 7), PectinateTree(7))),
    MASTSize(list(bal = BalancedTree(7), pec = PectinateTree(7)),
             as.phylo(0:3, 7)))
  expect_equal(t(NNIDist(BalancedTree(7), as.phylo(0:3, 7))),
    NNIDist(list(bal = BalancedTree(7), pec = PectinateTree(7)),
             as.phylo(0:3, 7))[1, , ])
})
