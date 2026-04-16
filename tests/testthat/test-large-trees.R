# Tests for large-tree support (> 2048 tips).
#
# These tests are guarded by skip_if(.SL_MAX_TIPS < required_tips) so they
# run only after TreeTools raises SL_MAX_TIPS and TreeDist is rebuilt.

test_that("R-level guard rejects trees exceeding .SL_MAX_TIPS", {
  skip_on_cran()
  too_many <- .SL_MAX_TIPS + 1L
  t1 <- TreeTools::RandomTree(too_many)
  t2 <- TreeTools::RandomTree(too_many)

  expect_error(ClusteringInfoDistance(t1, t2), "not yet supported")
  expect_error(PhylogeneticInfoDistance(t1, t2), "not yet supported")
  expect_error(MatchingSplitDistance(t1, t2), "not yet supported")
  expect_error(MatchingSplitInfoDistance(t1, t2), "not yet supported")
  expect_error(InfoRobinsonFoulds(t1, t2), "not yet supported")
  expect_error(NyeSimilarity(t1, t2), "not yet supported")
})

test_that("Batch path rejects trees exceeding .SL_MAX_TIPS", {
  skip_on_cran()
  too_many <- .SL_MAX_TIPS + 1L
  trees <- lapply(1:3, function(i) TreeTools::RandomTree(too_many))
  class(trees) <- "multiPhylo"

  expect_error(ClusteringInfoDistance(trees), "not yet supported")
  expect_error(PhylogeneticInfoDistance(trees), "not yet supported")
})

test_that("CID works for 4000-tip trees", {
  skip_on_cran()
  skip_if(.SL_MAX_TIPS < 4000L, "SL_MAX_TIPS not yet raised to 4000+")

  set.seed(7042)
  t1 <- TreeTools::RandomTree(4000)
  t2 <- TreeTools::RandomTree(4000)

  cid <- ClusteringInfoDistance(t1, t2)
  expect_type(cid, "double")
  expect_true(is.finite(cid))
  expect_gte(cid, 0)
})

test_that("Multiple metrics agree on identical 4000-tip trees", {
  skip_on_cran()
  skip_if(.SL_MAX_TIPS < 4000L, "SL_MAX_TIPS not yet raised to 4000+")

  set.seed(3891)
  t1 <- TreeTools::RandomTree(4000)

  expect_equal(ClusteringInfoDistance(t1, t1), 0)
  expect_equal(MatchingSplitDistance(t1, t1), 0)
  expect_equal(RobinsonFoulds(t1, t1), 0)
  expect_equal(InfoRobinsonFoulds(t1, t1), 0)
})

test_that("Batch CID works for 4000-tip trees", {
  skip_on_cran()
  skip_if(.SL_MAX_TIPS < 4000L, "SL_MAX_TIPS not yet raised to 4000+")

  set.seed(5283)
  trees <- lapply(1:5, function(i) TreeTools::RandomTree(4000))
  class(trees) <- "multiPhylo"

  d <- ClusteringInfoDistance(trees)
  expect_s3_class(d, "dist")
  expect_equal(attr(d, "Size"), 5L)
  expect_true(all(is.finite(d)))
  expect_true(all(d >= 0))
})

test_that("RF and IRF work for 8000-tip trees", {
  skip_on_cran()
  skip_if(.SL_MAX_TIPS < 8000L, "SL_MAX_TIPS not yet raised to 8000+")

  set.seed(6174)
  t1 <- TreeTools::RandomTree(8000)
  t2 <- TreeTools::RandomTree(8000)

  rf <- RobinsonFoulds(t1, t2)
  expect_type(rf, "double")
  expect_true(is.finite(rf))
  expect_gte(rf, 0)
  expect_equal(RobinsonFoulds(t1, t1), 0)

  irf <- InfoRobinsonFoulds(t1, t2)
  expect_type(irf, "double")
  expect_true(is.finite(irf))
  expect_gte(irf, 0)
})
