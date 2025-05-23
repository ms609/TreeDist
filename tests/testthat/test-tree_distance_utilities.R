library("TreeTools", quietly = TRUE)


test_that("Tree normalization works", {
  expect_equal(NormalizeInfo(unnormalized = 5,
                             tree1 = BalancedTree(3),
                             tree2 = BalancedTree(5),
                             how = TRUE,
                             InfoInTree = function(x, y) NTip(x) + NTip(y),
                             y = SingleTaxonTree()),
               5 / ((3 + 1) + (3 + 1)))
  expect_equal(NormalizeInfo(unnormalized = 0:4,
                             tree1 = lapply(1:5, BalancedTree),
                             tree2 = BalancedTree(3),
                             InfoInTree = NTip,
                             Combine = min),
               0:4 / c(1, 2, 3, 3, 3))
  expect_equal(NormalizeInfo(unnormalized = 3, 1, 1, how = 4L), 3 / 4)
  expect_equal(as.dist(matrix(0:24 / 10, 5, 5)),
               NormalizeInfo(as.dist(matrix(0:24, 5, 5)),
                             rep(5, 5), rep(5, 5), InfoInTree = I),
               ignore_attr = TRUE)
  expect_equal(matrix(1:20 / 10, 5, 4) * c(10,10,1,1,1, 10,10,1,1,1, rep(1,10)),
               NormalizeInfo(1:20, c(1, 1, 10, 10, 10), c(1, 1, 10, 10), 
                             InfoInTree = I, Combine = pmax))
})

test_that("CalculateTreeDistance() errs appropriately", {
  expect_error(CalculateTreeDistance(RobinsonFouldsSplits, "Not a tree"))
  expect_error(CalculateTreeDistance(RobinsonFouldsSplits, "Not a tree", BalancedTree(8)))
  expect_error(CalculateTreeDistance(RobinsonFouldsSplits, BalancedTree(8), "Not a tree"))
})

test_that("CalculateTreeDistance() handles splits appropriately", {
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
  expect_equal(dist(0), CalculateTreeDistance(RobinsonFouldsSplits, splits10, NULL))
  expect_equal(seq_len(7),
               attr(CalculateTreeDistance(RobinsonFouldsSplits, 
                                          splits10, splits10, TRUE), "matching")
               )
  expect_equal(CalculateTreeDistance(RobinsonFouldsSplits, splits10, splits10.1),
               CalculateTreeDistance(RobinsonFouldsSplits, splits10.1, splits10))
  expect_equal(CalculateTreeDistance(RobinsonFouldsSplits, splits10, splits10.3),
               CalculateTreeDistance(RobinsonFouldsSplits, splits10.3, splits10))
  
  expect_equal(ignore_attr = TRUE,
    CalculateTreeDistance(RobinsonFouldsSplits, splits10.3, splits10.3),
    CalculateTreeDistance(RobinsonFouldsSplits, splits10.3, trees10.3))
  
  expect_equal(ignore_attr = TRUE,
    as.matrix(CalculateTreeDistance(RobinsonFouldsSplits, trees10.3)),
    as.matrix(CalculateTreeDistance(RobinsonFouldsSplits, trees10.3, trees10.3)))
  
  expect_equal(ignore_attr = TRUE,
    as.matrix(CalculateTreeDistance(RobinsonFouldsSplits, trees10.3, trees10.3))[, c(1, 3, 2)],
    CalculateTreeDistance(RobinsonFouldsSplits, trees10.3, trees10.3[c(1, 3, 2)]))
  
  mat <- matrix(NA, 3, 4)
  for (i in 1:3) for (j in 1:4) {
    mat[i, j] <- CalculateTreeDistance(RobinsonFouldsSplits, splits10.3[[i]], splits10.4[[j]])
  }
  
  expect_equal(mat,
               CalculateTreeDistance(RobinsonFouldsSplits, splits10.3, splits10.4),
               ignore_attr = TRUE)
  
  expect_equal(ignore_attr = TRUE,
    CalculateTreeDistance(RobinsonFouldsSplits, splits10.4, splits10.3),
    t(CalculateTreeDistance(RobinsonFouldsSplits, splits10.3, splits10.4)))
})

test_that("Matches are reported", {
  # Trees copied from test-tree_distance.R
  treeSym8 <- ape::read.tree(text="((e, (f, (g, h))), (((a, b), c), d));")
  treeBal8 <- ape::read.tree(text="(((e, f), (g, h)), ((a, b), (c, d)));")
  treeTwoSplits <- ape::read.tree(text="(((a, b), c, d), (e, f, g, h));")
  
  A <- TRUE
  B <- FALSE
  splitsA <- as.Splits(c(A, A, B, B, B, B, B, B), tipLabels = letters[1:8])
  splitsB <- as.Splits(matrix(c(A, A, A, A, A, A, B, B, 
                          A, A, A, B, B, B, B, B,
                          B, B, A, A, A, A, A, A,
                          A, A, A, A, A, B, A, B), nrow = 4, byrow = TRUE),
                       tipLabels = letters[1:8])
  match <- TreeTools::match # Avoid being re-masked by "base"
  
  expect_equal("a b | c d e f g h => a b c | d e f g h",
               ReportMatching(splitsA, splitsB[[2]]))
  
  splits1 <- as.Splits(treeSym8)
  splits2 <- as.Splits(treeBal8)
  
  matchedSplits <- match(splits1, splits2)
  cs <- CompatibleSplits(splits1, splits2)
  cs[, matchedSplits] <- FALSE
  unmatched <- is.na(matchedSplits)
  matchedSplits[unmatched] <- apply(cs[unmatched, ], 1, which)
  expect_equal(matchedSplits,
               attr(GeneralizedRF(as.Splits(treeSym8), as.Splits(treeBal8), 8L, 
                                  cpp_shared_phylo, maximize = TRUE, 
                                  reportMatching = TRUE),
                    "matching"))
    
  expect_equal(matchedSplits,
               attr(SharedPhylogeneticInfoSplits(as.Splits(treeSym8),
                                                 as.Splits(treeBal8),
                                                 reportMatching = TRUE),
                    "matching"))
  
  tree1 <- PectinateTree(letters[1:8])
  tree2 <- BalancedTree(letters[8:1])
  splits1 <- as.Splits(tree1)
  splits2 <- as.Splits(tree2, tree1)
  
  .ExpectAtOK <- function(at) {
    expect_equal(4L, length(at))
    expect_equal(names(at),
                 c("matching", "matchedSplits", "matchedScores", "pairScores"))
  }
  
  Test <- function(Func, relaxed = FALSE, ...) {
    
    at <- attributes(Func(tree1, tree2, reportMatching = TRUE, ...))
    .ExpectAtOK(at)
    
    matchedSplits <- match(splits1, splits2)
    if (relaxed) {
      expect_equal(matchedSplits[!is.na(matchedSplits)],
                   as.integer(at[["matching"]][c(1, 3, 5)]))
    } else {
      cs <- CompatibleSplits(splits1, splits2)
      cs[, matchedSplits] <- FALSE
      unmatched <- is.na(matchedSplits)
      matchedSplits[unmatched] <- apply(cs[unmatched, ], 1, which)
      
      expect_equal(matchedSplits, as.integer(at$matching))
    }
    ghSplit <- at[["matchedSplits"]][
      match(as.Splits(c(rep(FALSE, 6), TRUE, TRUE), letters[1:8]),
            splits1[[which(!is.na(matchedSplits))]])]
    expect_equal(ghSplit, "g h | a b c d e f => g h | a b c d e f")
    
    at <- attributes(Func(treeSym8, treeTwoSplits, reportMatching = TRUE, ...))
    .ExpectAtOK(at)
    expect_equal(as.integer(at[["matching"]]),
                 match(as.Splits(treeSym8), as.Splits(treeTwoSplits, treeSym8)))
    expect_equal(at[["matchedSplits"]][[2]],
                 "a b | e f g h c d => a b | e f g h c d")
  }
  
  Test(SharedPhylogeneticInfo)
  Test(DifferentPhylogeneticInfo)
  Test(MatchingSplitInfo)
  Test(MatchingSplitInfoDistance)
  Test(MutualClusteringInfo)
  Test(ClusteringInfoDistance)
  
  Test(NyeSimilarity)
  Test(MatchingSplitDistance)
  Test(JaccardRobinsonFoulds, k = 2, allowConflict = FALSE)
  Test(JaccardRobinsonFoulds, k = 2, allowConflict = TRUE)
  Test(RobinsonFoulds, relaxed = TRUE)
  Test(InfoRobinsonFoulds, relaxed = TRUE)
  
  # Matching Split Distance matches differently:  
  at <- attributes(MatchingSplitDistance(treeSym8, treeBal8, 
                                         reportMatching = TRUE))
  .ExpectAtOK(at)
  expect_equal(as.integer(at$matching), c(1:3, 5:4))
  expect_equal(at[["matchedSplits"]][[5]],
               "a b | e f g h c d => a b | e f g h c d")
  
  # Zero match:
  expect_true(attr(SharedPhylogeneticInfo(
    ape::read.tree(text = "((a, b), (c, d));"),
    ape::read.tree(text = "((a, c), (b, d));"), 
    reportMatching = TRUE),
    "matchedSplits") %in% c(
      "a b | c d .. a c | b d",
      "a b | c d .. b d | a c",
      "c d | a b .. a c | b d",
      "c d | a b .. b d | a c"))
})

test_that("Matchings are calculated in both directions", {
  tree1 <- ape::read.tree(text="((1, 2), ((3, 4, (5, 9)), (6, (7, 8))));")
  tree2 <- ape::read.tree(text="((1, 2), ((3, (4, 5)), (6, (7, (8, 9)))));")
  splits1 <- as.Splits(tree1)
  splits2 <- as.Splits(tree2)
  nMatches <- min(length(splits1), length(splits2))
  
  Test <- function(Func, ...) {
    matching12 <- Func(tree1, tree2, reportMatching = TRUE, ...)
    expect_equal(nMatches, length(attr(matching12, "matchedSplits")))
    
    matching21 <- Func(tree2, tree1, reportMatching = TRUE, ...)
    expect_equal(nMatches, length(attr(matching21, "matchedSplits")))
  }
  
  Test(SharedPhylogeneticInfo)
  Test(DifferentPhylogeneticInfo)
  Test(MatchingSplitInfo)
  Test(MatchingSplitInfoDistance)
  Test(MutualClusteringInfo)
  Test(ClusteringInfoDistance)
  
  Test(NyeSimilarity)
  Test(JaccardRobinsonFoulds, k = 2, allowConflict = FALSE)
  Test(JaccardRobinsonFoulds, k = 2, allowConflict = TRUE)
  Test(MatchingSplitDistance)
  
  nMatches <- 1L
  Test(RobinsonFoulds)
  Test(InfoRobinsonFoulds)
  
})

test_that(".TreeDistance() supports all sizes", {
  expect_equal(rbind(
    bal = MASTSize(BalancedTree(7), as.phylo(0:3, 7)),
    pec = MASTSize(as.phylo(0:3, 7), PectinateTree(7))),
    MASTSize(list(bal = BalancedTree(7), pec = PectinateTree(7)),
             as.phylo(0:3, 7)))
  expect_equal(t(NNIDist(BalancedTree(7), as.phylo(0:3, 7))),
               NNIDist(list(bal = BalancedTree(7), pec = PectinateTree(7)),
                       as.phylo(0:3, 7))[1, , ])
  expect_error(.TreeDistance(RobinsonFoulds, PectinateTree(1:6),
                             PectinateTree(6)))
})

test_that("Unrooteds are handled by MAST", {
  trees <- list(UnrootTree(BalancedTree(8)), UnrootTree(PectinateTree(8)))
  expect_equal(6L, as.integer(MASTSize(trees, rooted = FALSE)),
               ignore_attr = TRUE)
})

test_that("Entropy() supports dots input", {
  expect_identical(2, Entropy(rep(1/4, 4)))
  expect_identical(1, Entropy(0, .5, 1/2))
})

test_that(".SharedOnly() works", {
  ah <- letters[1:8]
  bi <- letters[2:9]
  bh <- letters[2:8]
  ch <- letters[3:8]
  ci <- letters[3:9]
  balAH <- BalancedTree(ah)
  pecAH <- PectinateTree(ah)
  pecBI <- PectinateTree(bi)
  balCH <- BalancedTree(ch)
  pecCI <- PectinateTree(ci)
  
  expect_equal( # 1 v 1
    .SharedOnly(balAH, balCH),
    list(KeepTip(balAH, intersect(ah, ch)),
         KeepTip(balCH, intersect(ah, ch)))
  )
  expect_equal( # 1 v N
    .SharedOnly(balAH, c(balAH, balCH)),
    list(list(balAH, KeepTip(balAH, ch)),
         list(balAH, balCH))
  )
  expect_equal( # N v 1
    .SharedOnly(c(balAH, balCH), balAH),
    list(list(balAH, balCH),
         list(balAH, KeepTip(balAH, ch)))
  )
  expect_equal( # N v N, all tips differ
    .SharedOnly(c(balAH, balCH), c(pecBI, pecCI)),
    list(
      list(KeepTip(balAH, bh), KeepTip(balAH, ch),
           balCH, balCH),
      list(KeepTip(pecBI, bh), KeepTip(pecCI, ch),
           KeepTip(pecBI, ch), KeepTip(pecCI, ch))
    )
  )
  expect_equal( # N v N, each list same tips
    .SharedOnly(c(balAH, pecAH), c(pecBI, pecBI)),
    list(
      list(KeepTip(balAH, bh), KeepTip(balAH, bh),
           KeepTip(pecAH, bh), KeepTip(pecAH, bh)),
      list(KeepTip(pecBI, bh), KeepTip(pecBI, bh),
           KeepTip(pecBI, bh), KeepTip(pecBI, bh))
    )
  )
  expect_equal( # N v M, first tips same
    .SharedOnly(c(balAH, balAH), c(pecBI, pecCI, balAH)),
    list(
      list(KeepTip(balAH, bh), KeepTip(balAH, ch), balAH,
           KeepTip(balAH, bh), KeepTip(balAH, ch), balAH),
      list(KeepTip(pecBI, bh), KeepTip(pecCI, ch), balAH,
           KeepTip(pecBI, bh), KeepTip(pecCI, ch), balAH)
    )
  )
  expect_equal( # N v M, second tips same
    .SharedOnly(c(pecBI, pecCI, balAH), c(balAH, balAH)),
    list(
      list(KeepTip(pecBI, bh), KeepTip(pecBI, bh),
           KeepTip(pecCI, ch), KeepTip(pecCI, ch),
           balAH, balAH),
      list(KeepTip(balAH, bh), KeepTip(balAH, bh),
           KeepTip(balAH, ch), KeepTip(balAH, ch),
           balAH, balAH)
    )
  )
  expect_equal(.SharedOnly(balAH, NULL), list(NULL, NULL))
  expect_equal(
    .SharedOnly(c(balAH, pecBI, balCH), NULL),
    list(
      list(KeepTip(balAH, bh), KeepTip(balAH, ch),
           KeepTip(pecBI, ch)),
      list(KeepTip(pecBI, bh), balCH,
           balCH)
    )
  )
  
  t1 <- structure(list(
    t1_1 = structure(list(edge = structure(
      c(7L, 8L, 9L, 9L, 8L, 7L, 10L, 10L, 11L, 11L, 8L, 9L, 1L, 6L, 2L, 10L, 
        3L, 11L, 4L, 5L), dim = c(10L, 2L)),
      Nnode = 5L, tip.label = c("1", "2", "3", "4", "5", "6")),
      order = "preorder", class = "phylo"),
    t1_2 = structure(list(edge = structure(
      c(7L, 8L, 8L, 7L, 9L, 10L, 10L, 9L, 11L, 11L, 8L, 1L, 6L, 9L, 10L, 2L, 3L,
        11L, 4L, 5L), dim = c(10L, 2L)),
      Nnode = 5L, tip.label = c("1", "2", "3", "4", "5", "6")),
      order = "preorder", class = "phylo")),
    firstHit = c(seed = 0, final = 2), class = "multiPhylo")
  t2 <- structure(list(
    t2_1 = structure(list(
      edge = structure(c(6L, 6L, 7L, 7L, 8L, 8L, 9L, 9L, 1L, 7L, 2L, 8L, 3L, 9L,
                         4L, 5L), dim = c(8L, 2L)),
      Nnode = 4L, tip.label = c("2", "3", "4", "5", "6")), order = "preorder",
      class = "phylo")),
    firstHit = c(seed = 0, final = 1), class = "multiPhylo")
  expect_equal(class(.SharedOnly(t1, t2)[[2]][[1]]), class(t2[[1]]))
  expect_equal(class(.SharedOnly(t2, t1)[[1]][[1]]), class(t2[[1]]))
  expect_equal(class(.SharedOnly(t2, t2)[[1]][[1]]), class(t2[[1]]))
  
})
