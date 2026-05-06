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

test_that("ClusteringInfoDistance returns correct value for n > SL_MAX_TIPS", {
  skip_on_cran()
  n <- 2049L
  skip_if(TreeDist:::cpp_sl_max_tips() < n,
          "Requires TreeTools >= 2.3.0 (SL_MAX_TIPS > 2048)")

  t1 <- as.phylo(1, n)
  t2 <- as.phylo(2, n)

  # The individual-pair path and the all-pairs OpenMP batch path are
  # independent implementations; agreement confirms the correct value.
  expect_equal(
    ClusteringInfoDistance(t1, t2)[[1]],
    as.matrix(ClusteringInfoDistance(list(t1, t2)))[1, 2]
  )
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

test_that("Tip-count ceiling is enforced correctly", {
  skip_on_cran()
  skip_if(TreeDist:::cpp_sl_max_tips() < 2049L,
          "Requires TreeTools >= 2.3.0 (SL_MAX_TIPS > 2048)")

  expect_no_error(.CheckMaxTips(32767L))
  expect_error(.CheckMaxTips(32768L), "not yet supported")
})
