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
  
  expect_equal(c(5, 2, 3, 1, 4), attr(CGRF(as.Splits(treeSym8),
                                           as.Splits(treeBal8), 8L, 
                                           cpp_mutual_phylo,
                                           maximize = TRUE, 
                                           reportMatching = TRUE), 'matching'))
    
  expect_equal(c(5, 2, 3, 1, 4), attr(
    MutualPhylogeneticInfoSplits(as.Splits(treeSym8), as.Splits(treeBal8),
                               reportMatching = TRUE),
    'matching'))
  
  Test <- function (Func, relaxed = FALSE, ...) {
    
    at <- attributes(Func(PectinateTree(letters[1:8]),
                          BalancedTree(letters[8:1]), reportMatching = TRUE,
                          ...))
    expect_equal(3L, length(at))
    if (relaxed) {
      expect_equal(c(5L, 3L, 1L), as.integer(at$matching[c(1, 3, 5)]))
      ghSplit <- at$matchedSplits[3]
    } else {
      expect_equal(c(5L, 2L, 3L, 4L, 1L), as.integer(at$matching))
      ghSplit <- at$matchedSplits[5]
    }
    expect_equal('g h | a b c d e f => g h | a b c d e f', ghSplit)
    
    at <- attributes(Func(treeSym8, treeTwoSplits, reportMatching = TRUE,
                          ...))
    expect_equal(3L, length(at))
    expect_equal(c(NA, NA, 2L, NA, 1L), as.integer(at$matching))
    expect_equal('a b | e f g h c d => a b | e f g h c d', at$matchedSplits[2])
  }
  
  Test(MutualPhylogeneticInfo)
  Test(VariationOfPhylogeneticInfo)
  Test(MutualMatchingSplitInfo)
  Test(VariationOfMatchingSplitInfo)
  Test(MutualClusteringInfo)
  Test(VariationOfClusteringInfo)
  Test(MutualPhylogeneticInfo)
  Test(NyeTreeSimilarity)
  Test(MatchingSplitDistance)
  Test(JaccardRobinsonFoulds, k = 2, arboreal = FALSE)
  Test(JaccardRobinsonFoulds, k = 2, arboreal = TRUE)
  Test(RobinsonFoulds, relaxed = TRUE)
  Test(RobinsonFouldsInfo, relaxed = TRUE)

    # Matching Split Distance matches differently:  
  at <- attributes(MatchingSplitDistance(treeSym8, treeBal8, 
                                         reportMatching = TRUE))
  expect_equal(3L, length(at))
  expect_equal(c(1:3, 5:4), as.integer(at$matching))
  expect_equal('a b | e f g h c d => a b | e f g h c d', at$matchedSplits[5])
  
  # Zero match:
  expect_equal('c d | a b .. b d | a c', 
               attr(MutualPhylogeneticInfo( 
                 ape::read.tree(text="((a, b), (c, d));"),
                 ape::read.tree(text="((a, c), (b, d));"), 
                 reportMatching = TRUE), 'matchedSplits'))
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
  
  Test(MutualPhylogeneticInfo)
  Test(VariationOfPhylogeneticInfo)
  Test(MutualMatchingSplitInfo)
  Test(VariationOfMatchingSplitInfo)
  Test(MutualClusteringInfo)
  Test(VariationOfClusteringInfo)
  Test(MutualPhylogeneticInfo)
  Test(NyeTreeSimilarity)
  Test(JaccardRobinsonFoulds, k = 2, arboreal = FALSE)
  Test(JaccardRobinsonFoulds, k = 2, arboreal = TRUE)
  Test(MatchingSplitDistance)
  
  nMatches <- 1L
  Test(RobinsonFoulds)
  Test(RobinsonFouldsInfo)
  
})
