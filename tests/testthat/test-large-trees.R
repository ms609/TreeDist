test_that("Known-answer large-tree near-neighbours (4000 tips)", {
  skip_on_cran()

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
  expect_equal(unname(as.matrix(ClusteringInfoDistance(trees))[2, 1]),
               unname(cid), tolerance = 1e-8)
  expect_equal(unname(as.matrix(MatchingSplitDistance(trees))[2, 1]),
               unname(msd), tolerance = 1e-8)
  expect_equal(unname(as.matrix(InfoRobinsonFoulds(trees))[2, 1]),
               unname(irf), tolerance = 1e-8)
})
