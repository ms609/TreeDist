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
  splitsA <- matrix(c(A, A, B, B, B, B, B, B), ncol=1, 
                    dimnames = list(letters[1:8], NULL))
  splitsB <- matrix(rev(c(A, A, A, A, A, A, B, B, 
                          A, A, A, B, B, B, B, B,
                          B, B, A, A, A, A, A, A,
                          A, A, A, A, A, B, A, B)), ncol=4,
                    dimnames = list(letters[8:1], NULL))
  
  expect_equal("a b : c d e f g h => c d e f g h : a b",
               ReportMatching(splitsA, 
                              splitsB[rownames(splitsA), 2, drop=FALSE],
                              taxonNames = letters[1:8]))
  
  Test <- function (Func, relaxed = FALSE) {
    at <- attributes(Func(treeSym8, treeBal8, reportMatching = TRUE))
    expect_equal(3L, length(at))
    if (relaxed) {
      expect_equal(c(1L, 3L, 4L), as.integer(at$matching[c(1, 3, 5)]))
    } else {
      expect_equal(c(1L, 5L, 3L, 2L, 4L), as.integer(at$matching))
    }
    expect_equal('a b : e f g h c d => a b : e f g h c d', at$matchedSplits[5])
    
    at <- attributes(Func(treeSym8, treeTwoSplits, reportMatching = TRUE))
    expect_equal(3L, length(at))
    expect_equal(c(1L, NA, NA, NA, 2L), as.integer(at$matching))
    expect_equal('a b : e f g h c d => e f g h c d : a b', at$matchedSplits[2])
  }
  
  Test(MutualPhylogeneticInfo)
  Test(VariationOfPhylogeneticInfo)
  Test(MutualMatchingSplitInfo)
  Test(VariationOfMatchingSplitInfo)
  Test(MutualClusteringInfo)
  Test(VariationOfClusteringInfo)
  Test(MutualPhylogeneticInfo)
  Test(RobinsonFoulds, relaxed = TRUE)
  Test(NyeTreeSimilarity)

  # Matching Split Distance matches differently:  
  at <- attributes(MatchingSplitDistance(treeSym8, treeBal8, 
                                         reportMatching = TRUE))
  expect_equal(3L, length(at))
  expect_equal(c(1:3, 5:4), as.integer(at$matching))
  expect_equal('a b : e f g h c d => a b : e f g h c d', at$matchedSplits[5])
})
