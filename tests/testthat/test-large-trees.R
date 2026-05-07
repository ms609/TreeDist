library("TreeTools", quietly = TRUE)

# Tests that distance functions handle trees larger than the old SL_MAX_TIPS
# threshold (2048 in TreeTools <= 2.2.0).  These require TreeTools >= 2.3.0
# for heap-backed split storage; they skip gracefully on older builds.

test_that("Distance functions handle trees exceeding TT 2.2.0 limit", {
  skip_on_cran()
  n <- 2049L
  skip_if(TreeDist:::cpp_sl_max_tips() < n,
          "Requires TreeTools >= 2.3.0 (SL_MAX_TIPS > 2048)")

  t1 <- as.phylo(1, n)
  t2 <- as.phylo(2, n)

  expect_gt(RobinsonFoulds(t1, t2)[[1]], 0)
  expect_equal(RobinsonFoulds(t1, t1)[[1]], 0)

  expect_gt(InfoRobinsonFoulds(t1, t2)[[1]], 0)
  expect_equal(InfoRobinsonFoulds(t1, t1)[[1]], 0)

  expect_gt(ClusteringInfoDistance(t1, t2)[[1]], 0)
  expect_equal(ClusteringInfoDistance(t1, t1)[[1]], 0)

  expect_gt(DifferentPhylogeneticInfo(t1, t2)[[1]], 0)
  expect_equal(DifferentPhylogeneticInfo(t1, t1)[[1]], 0)

  expect_gt(MatchingSplitInfoDistance(t1, t2)[[1]], 0)
  expect_equal(MatchingSplitInfoDistance(t1, t1)[[1]], 0)
})

test_that("Distance scores agree across stack/heap storage threshold", {
  skip_on_cran()
  skip_if_not_installed("TreeTools", "2.3.0")
  # SL_STACK_SPLITS in TreeTools 2.3.0+ admits trees of n_tips <= 2048 to
  # stack-backed split storage; n_tips == 2049 forces the heap path.  This
  # test compares each metric on a stack-storage tree pair against the same
  # pair after inserting a cherry next to tip 1, which both pushes them
  # across the threshold and leaves their split structure invariant: each
  # existing split simply gains the new tip on whichever side already
  # contained tip 1, and one perfectly-matched cherry split is added
  # (contributing 0 to any distance).
  stackMax <- 2048L
  skip_if(TreeDist:::cpp_sl_max_tips() <= stackMax,
          "Requires TreeTools >= 2.3.0 (heap-backed storage)")

  t1s <- as.phylo(1234, stackMax)
  t2s <- PectinateTree(stackMax)
  t1h <- AddTip(t1s, where = 1, label = "extra")
  t2h <- AddTip(t2s, where = 1, label = "extra")

  # RF (count of unmatched splits) is invariant under the inserted cherry.
  expect_equal(RobinsonFoulds(t1s, t2s)[[1]],
               RobinsonFoulds(t1h, t2h)[[1]])

  # Per-split info depends on n_tips and on side sizes, so absolute
  # information-theoretic distances drift by a few % when n grows by 1;
  # normalized distances are far less sensitive.
  tol <- 0.002

  expect_equal(InfoRobinsonFoulds(t1s, t2s, normalize = TRUE)[[1]],
               InfoRobinsonFoulds(t1h, t2h, normalize = TRUE)[[1]],
               tolerance = tol)
  expect_equal(ClusteringInfoDistance(t1s, t2s, normalize = TRUE)[[1]],
               ClusteringInfoDistance(t1h, t2h, normalize = TRUE)[[1]],
               tolerance = tol)
  expect_equal(DifferentPhylogeneticInfo(t1s, t2s, normalize = TRUE)[[1]],
               DifferentPhylogeneticInfo(t1h, t2h, normalize = TRUE)[[1]],
               tolerance = tol)
  expect_equal(MatchingSplitInfoDistance(t1s, t2s, normalize = TRUE)[[1]],
               MatchingSplitInfoDistance(t1h, t2h, normalize = TRUE)[[1]],
               tolerance = tol)
})

test_that("RF and IRF work for 8000-tip trees", {
  skip_on_cran()
  skip_if(TreeDist:::cpp_sl_max_tips() < 2049L,
          "Requires TreeTools >= 2.3.0 (SL_MAX_TIPS > 2048)")
  skip_if(!getOption("slowMode", FALSE), "Only runs in slow mode")

  n <- 8000
  t1 <- PectinateTree(n)
  t2 <- BalancedTree(n)

  expect_gt(RobinsonFoulds(t1, t2)[[1]], 0)
  expect_equal(RobinsonFoulds(t1, t1)[[1]], 0)

  expect_gt(InfoRobinsonFoulds(t1, t2)[[1]], 0)
  expect_equal(InfoRobinsonFoulds(t1, t1)[[1]], 0)
})

test_that("Slow lookup path produces identical SPI to fast path", {
  skip_on_cran()
  # Coverage for the n_tips > SL_MAX_TIPS branches in shared_phylo:
  #   - lg2_*_slow fallback functions
  #   - one_overlap<false>, one_overlap_notb<false>, spi_overlap<false>
  # These are only triggered with > 32705 tips in production, but the
  # cpp_shared_phylo entry accepts a force_slow flag so the same code paths
  # can be exercised on small trees.  Both paths must produce identical
  # values: lg2_rooted_lookup() == lg2_rooted[] within the table range.
  t1 <- as.phylo(0, 30)
  t2 <- as.phylo(11, 30)
  splits1 <- as.Splits(t1)
  splits2 <- as.Splits(t2)

  fast <- TreeDist:::cpp_shared_phylo(splits1, splits2, 30L, force_slow = FALSE)
  slow <- TreeDist:::cpp_shared_phylo(splits1, splits2, 30L, force_slow = TRUE)

  expect_equal(fast$score, slow$score)
  expect_equal(fast$matching, slow$matching)

  # Identical-tree case: both paths must report the self-information score.
  fast_self <- TreeDist:::cpp_shared_phylo(splits1, splits1, 30L,
                                           force_slow = FALSE)
  slow_self <- TreeDist:::cpp_shared_phylo(splits1, splits1, 30L,
                                           force_slow = TRUE)
  expect_equal(fast_self$score, slow_self$score)
})

test_that("Tip-count ceiling is enforced correctly", {
  skip_on_cran()
  skip_if(TreeDist:::cpp_sl_max_tips() < 2049L,
          "Requires TreeTools >= 2.3.0 (SL_MAX_TIPS > 2048)")

  expect_no_error(.CheckMaxTips(32767L))
  expect_error(.CheckMaxTips(32768L), "not yet supported")
})
