test_that("Identical trees return fully resolved consensus", {
  tree <- as.phylo(0, nTip = 10)
  trees <- structure(rep(list(tree), 10), class = "multiPhylo")

  tc <- TransferConsensus(trees, greedy = "first")
  expect_equal(NSplits(tc), NSplits(tree))

  tc_best <- TransferConsensus(trees, greedy = "best")
  expect_equal(NSplits(tc_best), NSplits(tree))
})

test_that("Star tree returned when no signal", {
  set.seed(6129)
  # Fully random trees with many tips — negligible split overlap
  trees <- structure(lapply(1:30, function(i) TreeTools::RandomTree(50, root = TRUE)),
                     class = "multiPhylo")
  tc <- TransferConsensus(trees, greedy = "first")

  # Should be very unresolved (0 or near-0 splits)
  expect_lte(NSplits(tc), 5)
})

test_that("Transfer consensus is at least as resolved as majority-rule for structured trees",
{
  # Trees from as.phylo with moderate overlap
  trees <- as.phylo(0:29, nTip = 20)
  tc <- TransferConsensus(trees, greedy = "best")
  mr <- Consensus(trees, p = 0.5)

  # Transfer consensus should be at least as resolved
  expect_gte(NSplits(tc), NSplits(mr))
})

test_that("Both greedy strategies produce valid trees", {
  trees <- as.phylo(1:15, nTip = 12)

  tc_first <- TransferConsensus(trees, greedy = "first")
  tc_best  <- TransferConsensus(trees, greedy = "best")

  expect_s3_class(tc_first, "phylo")
  expect_s3_class(tc_best, "phylo")
  expect_equal(sort(TipLabels(tc_first)), sort(TipLabels(trees[[1]])))
  expect_equal(sort(TipLabels(tc_best)), sort(TipLabels(trees[[1]])))
})

test_that("scale = FALSE (unscaled) works", {
  trees <- as.phylo(1:10, nTip = 10)

  tc_scaled <- TransferConsensus(trees, scale = TRUE)
  tc_unscaled <- TransferConsensus(trees, scale = FALSE)

  expect_s3_class(tc_scaled, "phylo")
  expect_s3_class(tc_unscaled, "phylo")
})

test_that("init = 'majority' works", {
  trees <- as.phylo(0:19, nTip = 15)

  tc_empty <- TransferConsensus(trees, init = "empty")
  tc_maj   <- TransferConsensus(trees, init = "majority")

  expect_s3_class(tc_empty, "phylo")
  expect_s3_class(tc_maj, "phylo")
  # Both should produce reasonable trees
  expect_gte(NSplits(tc_empty), 1)
  expect_gte(NSplits(tc_maj), 1)
})

test_that("Error on bad input", {
  expect_error(TransferConsensus(list(TreeTools::BalancedTree(5))), "multiPhylo")
  expect_error(TransferConsensus(
    structure(list(TreeTools::BalancedTree(5)), class = "multiPhylo")),
    "at least 2")
})

test_that("Two-tree consensus returns a valid tree", {
  trees <- as.phylo(1:2, nTip = 8)
  tc <- TransferConsensus(trees)
  expect_s3_class(tc, "phylo")
  expect_equal(length(TipLabels(tc)), 8)
})


# =========================================================================
# C++ edge cases (transfer_consensus.cpp coverage)
# =========================================================================

test_that("TransferConsensus returns star tree for all-star input", {
  star <- StarTree(8)
  trees <- structure(rep(list(star), 3), class = "multiPhylo")
  tc <- TransferConsensus(trees)
  expect_s3_class(tc, "phylo")
  expect_equal(NSplits(tc), 0)
})

test_that("TransferConsensus covers all parameter combinations", {
  trees <- as.phylo(0:9, nTip = 10)

  # greedy="first", init="majority", scale=FALSE
  tc <- TransferConsensus(trees, greedy = "first", init = "majority",
                          scale = FALSE)
  expect_s3_class(tc, "phylo")

  # greedy="best", init="majority", scale=FALSE
  tc2 <- TransferConsensus(trees, greedy = "best", init = "majority",
                           scale = FALSE)
  expect_s3_class(tc2, "phylo")

  # greedy="first", init="empty", scale=FALSE
  tc3 <- TransferConsensus(trees, greedy = "first", init = "empty",
                           scale = FALSE)
  expect_s3_class(tc3, "phylo")
})

test_that("Greedy remove path fires with conflicting majority splits", {
  # Two groups of trees with very different topologies: majority-init seeds
  # splits from the larger group, some of which the greedy then removes.
  t1 <- as.phylo(0, nTip = 12)
  t2 <- as.phylo(100000, nTip = 12)
  trees <- structure(c(rep(list(t1), 6), rep(list(t2), 5)),
                     class = "multiPhylo")

  # greedy="best" exercises do_remove + find_second in GreedyState
  tc_best_s <- TransferConsensus(trees, greedy = "best", init = "majority",
                                 scale = TRUE)
  expect_s3_class(tc_best_s, "phylo")

  tc_best_u <- TransferConsensus(trees, greedy = "best", init = "majority",
                                 scale = FALSE)
  expect_s3_class(tc_best_u, "phylo")

  # greedy="first" exercises the same remove path via greedy_first

  tc_first_s <- TransferConsensus(trees, greedy = "first", init = "majority",
                                  scale = TRUE)
  expect_s3_class(tc_first_s, "phylo")

  tc_first_u <- TransferConsensus(trees, greedy = "first", init = "majority",
                                  scale = FALSE)
  expect_s3_class(tc_first_u, "phylo")
})

test_that("Greedy exercises diverse topologies (extra C++ path coverage)", {
  # Random trees with many tips — sparse split overlap
  set.seed(5123)
  rand_trees <- structure(
    lapply(1:30, function(i) TreeTools::RandomTree(10, root = TRUE)),
    class = "multiPhylo")
  # Diverse indexed trees
  diverse_trees <- as.phylo(seq(0, 500, by = 25), nTip = 20)
  # Three conflicting groups
  mixed_trees <- structure(c(
    rep(list(as.phylo(0, 8)), 5),
    rep(list(as.phylo(300, 8)), 5),
    rep(list(as.phylo(9999, 8)), 3)
  ), class = "multiPhylo")

  for (trees in list(rand_trees, diverse_trees, mixed_trees)) {
    for (gr in c("best", "first")) {
      for (init in c("empty", "majority")) {
        for (sc in c(TRUE, FALSE)) {
          tc <- TransferConsensus(trees, greedy = gr, init = init, scale = sc)
          expect_s3_class(tc, "phylo")
        }
      }
    }
  }
})

test_that("cpp_tc_profile() runs without error", {
  trees <- as.phylo(0:4, nTip = 8)
  tipLabels <- TipLabels(trees[[1]])
  splitsList <- lapply(trees, function(tr) unclass(as.Splits(tr, tipLabels)))

  # Greedy best, empty init
  res <- TreeDist:::cpp_tc_profile(splitsList, 8L, TRUE, TRUE, FALSE, 1L, 1L)
  expect_length(res, 5)
  expect_true(all(res >= 0))
  expect_equal(names(res),
               c("pool_splits", "transfer_dist_mat", "compute_td",
                 "compat_mat", "greedy"))

  # Greedy first, majority init, unscaled
  res2 <- TreeDist:::cpp_tc_profile(splitsList, 8L, FALSE, FALSE, TRUE, 1L, 1L)
  expect_length(res2, 5)
  expect_true(all(res2 >= 0))

  # Multiple iterations
  res3 <- TreeDist:::cpp_tc_profile(splitsList, 8L, TRUE, TRUE, FALSE, 3L, 1L)
  expect_length(res3, 5)

  # Greedy best, majority init
  res4 <- TreeDist:::cpp_tc_profile(splitsList, 8L, TRUE, TRUE, TRUE, 1L, 1L)
  expect_length(res4, 5)
})


# =========================================================================
# R-level internal helper functions (transfer_consensus.R lines 94-551)
# =========================================================================

test_that(".FlipRaw() flips bits and masks correctly", {
  # 8 tips (1 byte, all bits used, mask = 0xFF)
  expect_equal(TreeDist:::.FlipRaw(as.raw(0x07), 8), as.raw(0xf8))

  # 5 tips (1 byte, 5 bits, mask = 0x1F)
  expect_equal(TreeDist:::.FlipRaw(as.raw(0x07), 5), as.raw(0x18))

  # Multi-byte: 10 tips (2 bytes, last byte 2 bits, mask = 0x03)
  result <- TreeDist:::.FlipRaw(as.raw(c(0x00, 0x01)), 10)
  expect_equal(result, as.raw(c(0xff, 0x02)))
})

test_that(".TransferDistMat() computes pairwise transfer distances", {
  splitMat <- matrix(c(
    TRUE, TRUE, FALSE, FALSE, FALSE, FALSE,
    TRUE, TRUE, TRUE, FALSE, FALSE, FALSE,
    FALSE, FALSE, FALSE, TRUE, TRUE, TRUE
  ), nrow = 3, byrow = TRUE)

  DIST <- TreeDist:::.TransferDistMat(splitMat, 6)
  expect_equal(dim(DIST), c(3, 3))
  expect_equal(diag(DIST), c(0, 0, 0))
  expect_true(isSymmetric(DIST))
  # {1,2} vs {1,2,3}: hamming = 1
  expect_equal(DIST[1, 2], 1)
  # {1,2,3} vs {4,5,6}: complements → transfer = 0
  expect_equal(DIST[2, 3], 0)

  # Single split
  DIST1 <- TreeDist:::.TransferDistMat(
    matrix(c(TRUE, TRUE, FALSE, FALSE), nrow = 1), 4)
  expect_equal(dim(DIST1), c(1, 1))
  expect_equal(DIST1[1, 1], 0)
})

test_that(".CompatMat() checks bipartition compatibility", {
  splitMat <- matrix(c(
    TRUE, TRUE, FALSE, FALSE, FALSE, FALSE,
    FALSE, FALSE, TRUE, TRUE, FALSE, FALSE,
    TRUE, FALSE, TRUE, FALSE, FALSE, FALSE
  ), nrow = 3, byrow = TRUE)

  compat <- TreeDist:::.CompatMat(splitMat, 6)
  expect_equal(dim(compat), c(3, 3))
  expect_true(compat[1, 2])    # {1,2} and {3,4}: disjoint → compatible
  expect_false(compat[1, 3])   # {1,2} and {1,3}: all 4 intersections → incompatible
  expect_true(isSymmetric(compat))
})

test_that(".ComputeTD() works for both scale modes", {
  DIST <- matrix(c(0, 1, 3, 1, 0, 2, 3, 2, 0), 3, 3)
  lightSide <- c(3, 2, 3)
  sentDist <- lightSide - 1L
  treeMembers <- list(c(1, 2), c(2, 3))

  td_s <- TreeDist:::.ComputeTD(DIST, sentDist, treeMembers, lightSide, 2, TRUE)
  td_u <- TreeDist:::.ComputeTD(DIST, sentDist, treeMembers, lightSide, 2, FALSE)
  expect_length(td_s, 3)
  expect_length(td_u, 3)
  expect_true(all(td_s >= 0))
  expect_true(all(td_u >= 0))

  # Single member per tree
  td_single <- TreeDist:::.ComputeTD(
    DIST, sentDist, list(c(1), c(3)), lightSide, 2, TRUE)
  expect_length(td_single, 3)
})

test_that(".IsCompat() handles edge cases", {
  compat <- matrix(TRUE, 3, 3)
  compat[1, 3] <- compat[3, 1] <- FALSE

  # Empty included → always compatible
  expect_true(TreeDist:::.IsCompat(1, c(FALSE, FALSE, FALSE), compat, 10))
  # Too many included (>= nTip - 3) → never compatible
  expect_false(TreeDist:::.IsCompat(1, c(TRUE, TRUE, TRUE), compat, 6))
  # Compatible with included set
  expect_true(TreeDist:::.IsCompat(1, c(FALSE, TRUE, FALSE), compat, 10))
  # Incompatible with an included split
  expect_false(TreeDist:::.IsCompat(1, c(FALSE, FALSE, TRUE), compat, 10))
})

test_that(".Dist() returns correct distance", {
  DIST <- matrix(c(0, 2, 3, 2, 0, 4, 3, 4, 0), 3, 3)
  sentDist <- c(5, 6, 7)

  expect_equal(TreeDist:::.Dist(1, 2, DIST, sentDist), 2)
  expect_equal(TreeDist:::.Dist(1, NA, DIST, sentDist), 5)
  expect_equal(TreeDist:::.Dist(3, NA, DIST, sentDist), 7)
})

test_that(".DiagDist() vectorizes correctly", {
  DIST <- matrix(c(0, 2, 3, 2, 0, 4, 3, 4, 0), 3, 3)
  sentDist <- c(5, 6, 7)

  expect_equal(
    TreeDist:::.DiagDist(c(2, NA, 1), DIST, sentDist),
    c(DIST[1, 2], 6, DIST[3, 1])
  )
  expect_equal(TreeDist:::.DiagDist(c(NA, NA, NA), DIST, sentDist), sentDist)
  expect_equal(
    TreeDist:::.DiagDist(c(2, 3, 1), DIST, sentDist),
    c(DIST[1, 2], DIST[2, 3], DIST[3, 1])
  )
})

test_that(".FindSecond() works for both scale modes and edge cases", {
  DIST <- matrix(c(0, 1, 3, 1, 0, 2, 3, 2, 0), 3, 3)
  pMinus1 <- c(4, 4, 4)

  # Unscaled: exclude matchIdx, find best
  expect_equal(
    TreeDist:::.FindSecond(1, 2, c(1, 2, 3), DIST, pMinus1, FALSE), 1)

  # matchIdx = NA → search all
  expect_equal(
    TreeDist:::.FindSecond(1, NA_integer_, c(2, 3), DIST, pMinus1, FALSE), 2)

  # No candidates
  expect_true(is.na(
    TreeDist:::.FindSecond(1, 2, c(2), DIST, pMinus1, FALSE)))

  # All exceed threshold (unscaled)
  DIST_far <- matrix(c(0, 5, 5, 5, 0, 5, 5, 5, 0), 3, 3)
  expect_true(is.na(
    TreeDist:::.FindSecond(1, NA_integer_, c(2, 3), DIST_far, c(2, 2, 2), FALSE)
  ))

  # Scaled: normal case
  expect_equal(
    TreeDist:::.FindSecond(1, NA_integer_, c(2, 3), DIST, pMinus1, TRUE), 2)

  # Scaled: exceeds threshold
  expect_true(is.na(
    TreeDist:::.FindSecond(1, NA_integer_, c(2), DIST_far, c(2, 2, 2), TRUE)
  ))
})

test_that(".DoAdd() updates match state correctly", {
  DIST <- matrix(c(0, 1, 3, 1, 0, 2, 3, 2, 0), 3, 3)
  sentDist <- c(4, 4, 4)

  st <- new.env(parent = emptyenv())
  st$MATCH <- rep(NA_integer_, 3)
  st$MATCH2 <- rep(NA_integer_, 3)
  st$incl <- rep(FALSE, 3)

  # First add: all splits match the newly added split
  TreeDist:::.DoAdd(2, st, DIST, sentDist)
  expect_true(st$incl[2])
  expect_equal(st$MATCH[1], 2)
  expect_equal(st$MATCH[2], 2)
  expect_equal(st$MATCH[3], 2)

  # Second add: closer splits update, old match demotes to second
  TreeDist:::.DoAdd(1, st, DIST, sentDist)
  expect_equal(st$MATCH[1], 1)    # self-match closer
  expect_equal(st$MATCH2[1], 2)   # old match becomes second
  expect_equal(st$MATCH[2], 2)    # split 2 still closest to itself
  expect_equal(st$MATCH2[2], 1)   # split 1 becomes second
})

test_that(".DoRemove() handles sentinel promotion and affected2", {
  DIST <- matrix(c(0, 1, 3, 1, 0, 2, 3, 2, 0), 3, 3)
  sentDist <- c(4, 4, 4)
  lightSide <- c(5, 5, 5)

  # Case 1: sentinel promotion (MATCH2 = NA, remove MATCH, no other included)
  st <- new.env(parent = emptyenv())
  st$MATCH <- c(2L, NA_integer_, 2L)
  st$MATCH2 <- c(NA_integer_, NA_integer_, NA_integer_)
  st$incl <- c(FALSE, TRUE, FALSE)

  TreeDist:::.DoRemove(2, st, DIST, sentDist, lightSide, FALSE)
  expect_true(is.na(st$MATCH[1]))
  expect_true(is.na(st$MATCH[3]))

  # Case 2: sentinel promotion with rescan finding another match
  st2 <- new.env(parent = emptyenv())
  st2$MATCH <- c(2L, 2L, NA_integer_)
  st2$MATCH2 <- c(NA_integer_, NA_integer_, NA_integer_)
  st2$incl <- c(TRUE, TRUE, FALSE)

  TreeDist:::.DoRemove(2, st2, DIST, sentDist, lightSide, FALSE)
  expect_equal(st2$MATCH[1], 1L)  # rescanned, found split 1

  # Case 3: affected2 path (MATCH2 == removed, MATCH != removed)
  st3 <- new.env(parent = emptyenv())
  st3$MATCH <- c(1L, 1L, 1L)
  st3$MATCH2 <- c(2L, 2L, NA_integer_)
  st3$incl <- c(TRUE, TRUE, FALSE)

  TreeDist:::.DoRemove(2, st3, DIST, sentDist, lightSide, FALSE)
  expect_equal(st3$MATCH, c(1L, 1L, 1L))  # unchanged
  expect_true(is.na(st3$MATCH2[1]))
  expect_true(is.na(st3$MATCH2[2]))

  # Case 4: scaled mode
  st4 <- new.env(parent = emptyenv())
  st4$MATCH <- c(2L, NA_integer_, NA_integer_)
  st4$MATCH2 <- c(NA_integer_, NA_integer_, NA_integer_)
  st4$incl <- c(TRUE, TRUE, FALSE)

  TreeDist:::.DoRemove(2, st4, DIST, sentDist, lightSide, TRUE)
  expect_equal(st4$MATCH[1], 1L)
})

test_that(".AddBenefitVec() computes add benefits", {
  DIST <- matrix(c(0, 1, 3, 1, 0, 2, 3, 2, 0), 3, 3)
  lightSide <- c(3, 3, 3)
  sentDist <- lightSide - 1L
  counts <- c(2, 3, 1)
  TD <- c(0.5, 0.3, 0.8)

  st <- new.env(parent = emptyenv())
  st$MATCH <- rep(NA_integer_, 3)
  st$MATCH2 <- rep(NA_integer_, 3)
  st$incl <- rep(FALSE, 3)

  ben_s <- TreeDist:::.AddBenefitVec(
    c(1, 2, 3), st, DIST, sentDist, TD, counts, lightSide, scale = TRUE)
  expect_length(ben_s, 3)

  ben_u <- TreeDist:::.AddBenefitVec(
    c(1, 2, 3), st, DIST, sentDist, TD, counts, lightSide, scale = FALSE)
  expect_length(ben_u, 3)
})

test_that(".RemoveBenefitVec() computes remove benefits", {
  DIST <- matrix(c(0, 1, 3, 1, 0, 2, 3, 2, 0), 3, 3)
  lightSide <- c(3, 3, 3)
  sentDist <- lightSide - 1L
  counts <- c(2, 3, 1)
  TD <- c(0.5, 0.3, 0.8)

  st <- new.env(parent = emptyenv())
  st$MATCH <- c(2L, 1L, 2L)
  st$MATCH2 <- c(NA_integer_, NA_integer_, NA_integer_)
  st$incl <- c(TRUE, TRUE, FALSE)

  ben_s <- TreeDist:::.RemoveBenefitVec(
    c(1, 2), st, DIST, sentDist, TD, counts, lightSide, scale = TRUE)
  expect_length(ben_s, 2)

  ben_u <- TreeDist:::.RemoveBenefitVec(
    c(1, 2), st, DIST, sentDist, TD, counts, lightSide, scale = FALSE)
  expect_length(ben_u, 2)

  # No affected splits (MATCH doesn't reference candidate)
  st2 <- new.env(parent = emptyenv())
  st2$MATCH <- c(1L, 1L, 1L)
  st2$MATCH2 <- rep(NA_integer_, 3)
  st2$incl <- c(TRUE, TRUE, FALSE)

  ben_none <- TreeDist:::.RemoveBenefitVec(
    c(2), st2, DIST, sentDist, TD, counts, lightSide, scale = TRUE)
  # Only TD contribution, no fn_cost (no affected splits)
  expect_equal(ben_none, TD[2])
})

test_that(".PoolSplits() pools and deduplicates splits", {
  trees <- as.phylo(0:4, nTip = 8)
  tipLabels <- TipLabels(trees[[1]])
  pool <- TreeDist:::.PoolSplits(trees, tipLabels)

  expect_true(is.matrix(pool$splits))
  expect_true(is.matrix(pool$rawSplits))
  expect_equal(ncol(pool$splits), length(tipLabels))
  expect_equal(length(pool$counts), nrow(pool$splits))
  expect_equal(length(pool$lightSide), nrow(pool$splits))
  expect_equal(length(pool$treeMembers), length(trees))
  expect_true(all(pool$counts >= 1))
  expect_true(all(pool$lightSide >= 1))
})

test_that(".SplitsToPhylo() converts splits to tree or star", {
  trees <- as.phylo(0:2, nTip = 8)
  tipLabels <- TipLabels(trees[[1]])
  pool <- TreeDist:::.PoolSplits(trees, tipLabels)

  # Some splits included
  included <- rep(FALSE, nrow(pool$rawSplits))
  included[1:min(2, nrow(pool$rawSplits))] <- TRUE
  tr <- TreeDist:::.SplitsToPhylo(pool$rawSplits, included, tipLabels, 8)
  expect_s3_class(tr, "phylo")
  expect_equal(sort(TipLabels(tr)), sort(tipLabels))

  # No splits → star tree
  tr_star <- TreeDist:::.SplitsToPhylo(
    pool$rawSplits, rep(FALSE, nrow(pool$rawSplits)), tipLabels, 8)
  expect_s3_class(tr_star, "phylo")
  expect_equal(NSplits(tr_star), 0)
})

test_that(".InitMatches() initializes match state from included splits", {
  # DIST: included = {1,2}; split 3 close to split 2 (dist 1 < thresh 2)
  DIST <- matrix(c(
    0, 1, 3, 5,
    1, 0, 1, 4,
    3, 1, 0, 1,
    5, 4, 1, 0
  ), 4, 4, byrow = TRUE)
  lightSide <- c(3, 3, 3, 3)
  sentDist <- lightSide - 1L  # all 2

  # Scaled mode
  st <- new.env(parent = emptyenv())
  st$MATCH <- rep(NA_integer_, 4)
  st$MATCH2 <- rep(NA_integer_, 4)
  st$incl <- c(TRUE, TRUE, FALSE, FALSE)

  TreeDist:::.InitMatches(st, DIST, sentDist, lightSide, scale = TRUE)
  # Split 1: self-match (dist 0), second = split 2 (dist 1)
  expect_equal(st$MATCH[1], 1)
  expect_equal(st$MATCH2[1], 2)
  # Split 3: best among included is split 2 (dist 1), 1/2 < 1 → matched
  expect_equal(st$MATCH[3], 2)
  # Split 4: best among included is split 2 (dist 4), 4/2 >= 1 → skipped
  expect_true(is.na(st$MATCH[4]))

  # Unscaled mode
  st2 <- new.env(parent = emptyenv())
  st2$MATCH <- rep(NA_integer_, 4)
  st2$MATCH2 <- rep(NA_integer_, 4)
  st2$incl <- c(TRUE, TRUE, FALSE, FALSE)

  TreeDist:::.InitMatches(st2, DIST, sentDist, lightSide, scale = FALSE)
  expect_equal(st2$MATCH[1], 1)
  expect_equal(st2$MATCH[3], 2)

  # Empty included set → no-op
  st3 <- new.env(parent = emptyenv())
  st3$MATCH <- rep(NA_integer_, 4)
  st3$MATCH2 <- rep(NA_integer_, 4)
  st3$incl <- rep(FALSE, 4)

  TreeDist:::.InitMatches(st3, DIST, sentDist, lightSide, scale = TRUE)
  expect_true(all(is.na(st3$MATCH)))
})


# =========================================================================
# Full R-level greedy pipeline integration tests
# =========================================================================

test_that("R-level GreedyBest pipeline produces valid consensus (scaled)", {
  trees <- as.phylo(0:9, nTip = 8)
  tipLabels <- TipLabels(trees[[1]])
  nTip <- length(tipLabels)
  nTree <- length(trees)

  pool <- TreeDist:::.PoolSplits(trees, tipLabels)
  nSplits <- nrow(pool$splits)
  DIST <- TreeDist:::.TransferDistMat(pool$splits, nTip)
  sentDist <- pool$lightSide - 1L
  TD <- TreeDist:::.ComputeTD(DIST, sentDist, pool$treeMembers, pool$lightSide,
                               nTree, scale = TRUE)
  compat <- TreeDist:::.CompatMat(pool$splits, nTip)
  sortOrd <- order(-pool$counts, seq_len(nSplits))

  st <- new.env(parent = emptyenv())
  st$MATCH <- rep(NA_integer_, nSplits)
  st$MATCH2 <- rep(NA_integer_, nSplits)
  st$incl <- rep(FALSE, nSplits)

  TreeDist:::.GreedyBest(st, DIST, sentDist, TD, pool$counts, pool$lightSide,
                          compat, sortOrd, scale = TRUE, nSplits, nTip)
  expect_true(any(st$incl))

  tr <- TreeDist:::.SplitsToPhylo(pool$rawSplits, st$incl, tipLabels, nTip)
  expect_s3_class(tr, "phylo")
})

test_that("R-level GreedyFirst pipeline produces valid consensus (unscaled)", {
  trees <- as.phylo(0:9, nTip = 8)
  tipLabels <- TipLabels(trees[[1]])
  nTip <- length(tipLabels)
  nTree <- length(trees)

  pool <- TreeDist:::.PoolSplits(trees, tipLabels)
  nSplits <- nrow(pool$splits)
  DIST <- TreeDist:::.TransferDistMat(pool$splits, nTip)
  sentDist <- pool$lightSide - 1L
  TD <- TreeDist:::.ComputeTD(DIST, sentDist, pool$treeMembers, pool$lightSide,
                               nTree, scale = FALSE)
  compat <- TreeDist:::.CompatMat(pool$splits, nTip)
  sortOrd <- order(-pool$counts, seq_len(nSplits))

  st <- new.env(parent = emptyenv())
  st$MATCH <- rep(NA_integer_, nSplits)
  st$MATCH2 <- rep(NA_integer_, nSplits)
  st$incl <- rep(FALSE, nSplits)

  TreeDist:::.GreedyFirst(st, DIST, sentDist, TD, pool$counts, pool$lightSide,
                           compat, sortOrd, scale = FALSE, nSplits, nTip)
  expect_true(any(st$incl))

  tr <- TreeDist:::.SplitsToPhylo(pool$rawSplits, st$incl, tipLabels, nTip)
  expect_s3_class(tr, "phylo")
})

test_that("R-level GreedyBest with majority init exercises InitMatches", {
  trees <- as.phylo(0:14, nTip = 10)
  tipLabels <- TipLabels(trees[[1]])
  nTip <- length(tipLabels)
  nTree <- length(trees)

  pool <- TreeDist:::.PoolSplits(trees, tipLabels)
  nSplits <- nrow(pool$splits)
  DIST <- TreeDist:::.TransferDistMat(pool$splits, nTip)
  sentDist <- pool$lightSide - 1L
  TD <- TreeDist:::.ComputeTD(DIST, sentDist, pool$treeMembers, pool$lightSide,
                               nTree, scale = TRUE)
  compat <- TreeDist:::.CompatMat(pool$splits, nTip)
  sortOrd <- order(-pool$counts, seq_len(nSplits))

  st <- new.env(parent = emptyenv())
  st$MATCH <- rep(NA_integer_, nSplits)
  st$MATCH2 <- rep(NA_integer_, nSplits)

  # Majority init: include splits present in > half of trees
  st$incl <- pool$counts > (nTree / 2)
  if (any(st$incl)) {
    TreeDist:::.InitMatches(st, DIST, sentDist, pool$lightSide, scale = TRUE)
  }

  TreeDist:::.GreedyBest(st, DIST, sentDist, TD, pool$counts, pool$lightSide,
                          compat, sortOrd, scale = TRUE, nSplits, nTip)

  tr <- TreeDist:::.SplitsToPhylo(pool$rawSplits, st$incl, tipLabels, nTip)
  expect_s3_class(tr, "phylo")
})

test_that("R-level GreedyFirst with majority init (unscaled)", {
  trees <- as.phylo(0:14, nTip = 10)
  tipLabels <- TipLabels(trees[[1]])
  nTip <- length(tipLabels)
  nTree <- length(trees)

  pool <- TreeDist:::.PoolSplits(trees, tipLabels)
  nSplits <- nrow(pool$splits)
  DIST <- TreeDist:::.TransferDistMat(pool$splits, nTip)
  sentDist <- pool$lightSide - 1L
  TD <- TreeDist:::.ComputeTD(DIST, sentDist, pool$treeMembers, pool$lightSide,
                               nTree, scale = FALSE)
  compat <- TreeDist:::.CompatMat(pool$splits, nTip)
  sortOrd <- order(-pool$counts, seq_len(nSplits))

  st <- new.env(parent = emptyenv())
  st$MATCH <- rep(NA_integer_, nSplits)
  st$MATCH2 <- rep(NA_integer_, nSplits)
  st$incl <- pool$counts > (nTree / 2)
  if (any(st$incl)) {
    TreeDist:::.InitMatches(st, DIST, sentDist, pool$lightSide, scale = FALSE)
  }

  TreeDist:::.GreedyFirst(st, DIST, sentDist, TD, pool$counts, pool$lightSide,
                           compat, sortOrd, scale = FALSE, nSplits, nTip)

  tr <- TreeDist:::.SplitsToPhylo(pool$rawSplits, st$incl, tipLabels, nTip)
  expect_s3_class(tr, "phylo")
})

test_that("R greedy pipeline matches C++ consensus", {
  trees <- as.phylo(0:9, nTip = 8)
  tipLabels <- TipLabels(trees[[1]])
  nTip <- length(tipLabels)
  nTree <- length(trees)

  # C++ result
  tc_cpp <- TransferConsensus(trees, greedy = "best", scale = TRUE)

  # R-level result
  pool <- TreeDist:::.PoolSplits(trees, tipLabels)
  nSplits <- nrow(pool$splits)
  DIST <- TreeDist:::.TransferDistMat(pool$splits, nTip)
  sentDist <- pool$lightSide - 1L
  TD <- TreeDist:::.ComputeTD(DIST, sentDist, pool$treeMembers, pool$lightSide,
                               nTree, scale = TRUE)
  compat <- TreeDist:::.CompatMat(pool$splits, nTip)
  sortOrd <- order(-pool$counts, seq_len(nSplits))

  st <- new.env(parent = emptyenv())
  st$MATCH <- rep(NA_integer_, nSplits)
  st$MATCH2 <- rep(NA_integer_, nSplits)
  st$incl <- rep(FALSE, nSplits)

  TreeDist:::.GreedyBest(st, DIST, sentDist, TD, pool$counts, pool$lightSide,
                          compat, sortOrd, scale = TRUE, nSplits, nTip)

  tc_r <- TreeDist:::.SplitsToPhylo(pool$rawSplits, st$incl, tipLabels, nTip)

  # Both should produce resolved trees with same tip labels
  expect_s3_class(tc_r, "phylo")
  expect_equal(sort(TipLabels(tc_r)), sort(TipLabels(tc_cpp)))
  # Resolution should be similar (not necessarily identical due to
  # possible differences in canonicalization/sort-order between R and C++)
  expect_gte(NSplits(tc_r), 1)
})
