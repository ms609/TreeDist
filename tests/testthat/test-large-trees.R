# Tests for large-tree support (> 2048 tips).
#
# These tests are guarded by skip_if(.SL_MAX_TIPS < required_tips) so they
# run only after TreeTools raises SL_MAX_TIPS and TreeDist is rebuilt.

test_that("R-level guard rejects trees exceeding .SL_MAX_TIPS", {
  skip_on_cran()
  too_many <- .SL_MAX_TIPS + 1L
  t1 <- as.phylo(0, too_many)
  t2 <- as.phylo(1, too_many)

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
  trees <- as.phylo(0:2, too_many)
  class(trees) <- "multiPhylo"

  expect_error(ClusteringInfoDistance(trees), "not yet supported")
  expect_error(PhylogeneticInfoDistance(trees), "not yet supported")
})

test_that("Known-answer large-tree near-neighbours (4000 tips)", {
  skip_on_cran()
  skip_if(.SL_MAX_TIPS < 4000L, "SL_MAX_TIPS not yet raised to 4000+")

  # Similar deterministic trees exercise shortcut paths and run quickly.
  t1 <- as.phylo(0, 4000)
  t2 <- as.phylo(1, 4000)
  trees <- structure(list(t1, t2), class = "multiPhylo")

  # Known answer for adjacent `as.phylo()` trees: one non-shared split per tree.
  rf <- RobinsonFoulds(t1, t2)
  expect_equal(rf, 2)
  expect_equal(as.matrix(RobinsonFoulds(trees))[2, 1], 2)

  # Other large-tree metrics should be finite and non-negative.
  cid <- ClusteringInfoDistance(t1, t2)
  msd <- MatchingSplitDistance(t1, t2)
  irf <- InfoRobinsonFoulds(t1, t2)
  expect_true(all(is.finite(c(cid, msd, irf))))
  expect_true(all(c(cid, msd, irf) >= 0))

  # Batch and pairwise paths must agree.
  expect_equal(as.matrix(ClusteringInfoDistance(trees))[2, 1], cid, tolerance = 1e-10)
  expect_equal(as.matrix(MatchingSplitDistance(trees))[2, 1], msd, tolerance = 1e-10)
  expect_equal(as.matrix(InfoRobinsonFoulds(trees))[2, 1], irf, tolerance = 1e-10)
})
