test_that("Known-answer large-tree near-neighbours (>2048 tips)", {

  # Similar deterministic trees exercise shortcut paths and run quickly.
  t1 <- as.phylo(0, 2050)
  t2 <- as.phylo(1, 2050)
  trees <- structure(list(t1, t2), class = "multiPhylo")

  # Known answer for adjacent `as.phylo()` trees: one non-shared split per tree.
  rf <- RobinsonFoulds(t1, t2)
  expect_equal(rf, 2)
  expect_equal(as.matrix(RobinsonFoulds(trees))[2, 1], 2)

  # Other large-tree metrics should be finite and non-negative.
  cid <- ClusteringInfoDistance(t1, t2)
  expect_equal(cid, 0.01409, tolerance = 1e-5)
  msd <- MatchingSplitDistance(t1, t2)
  expect_equal(msd, 2)
  irf <- InfoRobinsonFoulds(t1, t2)
  expect_equal(irf, 23.999, tolerance = 1e-4)
  
  # Batch and pairwise paths must agree.
  expect_equal(unname(as.matrix(ClusteringInfoDistance(trees))[2, 1]),
               unname(cid), tolerance = 1e-8)
  expect_equal(unname(as.matrix(MatchingSplitDistance(trees))[2, 1]),
               unname(msd), tolerance = 1e-8)
  expect_equal(unname(as.matrix(InfoRobinsonFoulds(trees))[2, 1]),
               unname(irf), tolerance = 1e-8)
})

test_that("Large-tree (>SL_MAX_TIPS) non-table paths: PID, MSID, Jaccard, MCI", {
  # 2052 tips: SL_MAX_TIPS (2048) is exceeded, forcing fallback (non-table)
  # scoring paths in shared_phylo, msi_distance and the batch equivalents.
  # Split size can reach 2051 > SL_MAX_TIPS + 1, triggering lg2_rooted_lookup
  # fallback.  Near-neighbour trees share most splits so the LAP is tiny.
  t1 <- as.phylo(0, 2052)
  t2 <- as.phylo(1, 2052)
  trees <- structure(list(t1, t2), class = "multiPhylo")

  pid <- PhylogeneticInfoDistance(t1, t2)
  expect_type(pid, "double")
  expect_true(is.finite(pid))
  expect_gt(pid, 0)

  msid <- MatchingSplitInfoDistance(t1, t2)
  expect_type(msid, "double")
  expect_true(is.finite(msid))
  expect_gt(msid, 0)

  # Batch and pairwise paths must agree for PID and MSID
  expect_equal(unname(as.matrix(PhylogeneticInfoDistance(trees))[2, 1]),
               unname(pid), tolerance = 1e-8)
  expect_equal(unname(as.matrix(MatchingSplitInfoDistance(trees))[2, 1]),
               unname(msid), tolerance = 1e-8)

  # NyeSimilarity exercises jaccard_similarity serial interrupt path
  jac <- NyeSimilarity(t1, t2)
  expect_type(jac, "double")
  expect_true(is.finite(jac))
  expect_gt(jac, 0)

  skip_on_cran()
  skip_if_not(getOption("slowMode", FALSE))
  # reportMatching = TRUE forces the serial mutual_clustering path with interrupt
  cid_match <- ClusteringInfoDistance(t1, t2, reportMatching = TRUE)
  expect_true(is.integer(attr(cid_match, "matching")))
  expect_gt(length(attr(cid_match, "matching")), 0)
})
