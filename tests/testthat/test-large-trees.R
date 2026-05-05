library("TreeTools", quietly = TRUE)

# Tests that distance functions handle trees larger than the old SL_MAX_TIPS
# threshold (2048 in TreeTools <= 2.2.0).  These require TreeTools >= 2.3.0
# for heap-backed split storage; they skip gracefully on older builds.

test_that("RF and IRF work for 4000-tip trees", {
  skip_on_cran()
  skip_if(TreeDist:::.SL_MAX_TIPS < 4000,
          "Requires TreeTools >= 2.3.0 (SL_MAX_TIPS >= 4000)")

  n <- 4000
  t1 <- PectinateTree(n)
  t2 <- BalancedTree(n)

  rf <- RobinsonFoulds(t1, t2)[[1]]
  expect_gt(rf, 0)
  expect_equal(RobinsonFoulds(t1, t1)[[1]], 0)

  irf <- InfoRobinsonFoulds(t1, t2)[[1]]
  expect_gt(irf, 0)
  expect_equal(InfoRobinsonFoulds(t1, t1)[[1]], 0)
})

test_that("RF and IRF work for 8000-tip trees", {
  skip_on_cran()
  skip_if(TreeDist:::.SL_MAX_TIPS < 4000,
          "Requires TreeTools >= 2.3.0 (SL_MAX_TIPS >= 4000)")
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
  skip_if(TreeDist:::.SL_MAX_TIPS < 4000,
          "Requires TreeTools >= 2.3.0")

  expect_no_error(.CheckMaxTips(32767L))
  expect_error(.CheckMaxTips(32768L), "not yet supported")
})
