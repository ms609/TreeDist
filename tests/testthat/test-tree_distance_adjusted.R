# Tests for AdjustedClusteringInfoDistance / AdjustedPhylogeneticInfoDistance
# and the underlying ExpectedClusteringInfoDistance /
# ExpectedPhylogeneticInfoDistance back-ends (lookup / enumerate / mc).
#
# The shipped lookup at n = 4..7 stores values from exact enumeration
# (`data-raw/enumerate_small_n.R`), so the lookup and enumerate code paths
# should agree to numerical precision there; the MC fallback at larger n
# should agree to within a few standard errors.

suppressPackageStartupMessages({
  library("TreeTools", quietly = TRUE)
})

test_that("ACID/APhID return 0 for identical trees", {
  for (n in c(6L, 20L, 50L)) {
    trees <- list(
      balanced = BalancedTree(n),
      pectinate = PectinateTree(n)
    )
    set.seed(n)
    trees[["random"]] <- RandomTree(n, root = TRUE)
    for (tr in trees) {
      expect_equal(AdjustedClusteringInfoDistance(tr, tr), 0)
      expect_equal(AdjustedPhylogeneticInfoDistance(tr, tr), 0)
    }
  }
})

test_that("ACID/APhID are symmetric", {
  n <- 20L
  set.seed(1)
  t1 <- RandomTree(n, root = TRUE)
  t2 <- RandomTree(n, root = TRUE)
  expect_equal(AdjustedClusteringInfoDistance(t1, t2),
               AdjustedClusteringInfoDistance(t2, t1))
  expect_equal(AdjustedPhylogeneticInfoDistance(t1, t2),
               AdjustedPhylogeneticInfoDistance(t2, t1))
})

test_that("Lookup and enumerate agree at n = 6 and n = 7", {
  # After splicing exact values into the lookup at n = 4..7, the two code
  # paths read the same numbers; this is a self-consistency check that both
  # paths are wired up correctly.
  for (n in c(6L, 7L)) {
    expect_equal(ExpectedClusteringInfoDistance(n),
                 ExpectedClusteringInfoDistance(n, method = "enumerate"),
                 tolerance = 1e-10)
    expect_equal(ExpectedPhylogeneticInfoDistance(n),
                 ExpectedPhylogeneticInfoDistance(n, method = "enumerate"),
                 tolerance = 1e-10)
  }
})

test_that("Lookup agrees with on-the-fly MC at n = 20", {
  n <- 20L
  nSim <- 2000L
  set.seed(20)
  mc_cid <- ExpectedClusteringInfoDistance(n, method = "mc", nSim = nSim)
  set.seed(20)
  mc_pid <- ExpectedPhylogeneticInfoDistance(n, method = "mc", nSim = nSim)
  ref_cid <- ExpectedClusteringInfoDistance(n)
  ref_pid <- ExpectedPhylogeneticInfoDistance(n)
  ref <- randomTreeDistances[randomTreeDistances[["n"]] == n, , drop = FALSE]
  # 3 SE tolerance using the stored `sd`.
  tol_cid <- 3 * ref[["cid_sd"]] / sqrt(nSim)
  tol_pid <- 3 * ref[["pid_sd"]] / sqrt(nSim)
  expect_equal(mc_cid, ref_cid, tolerance = tol_cid)
  expect_equal(mc_pid, ref_pid, tolerance = tol_pid)
})

test_that("Lookup falls back to MC outside the shipped range", {
  # n = 500 is beyond the shipped table's max (n = 400 post-extension); the
  # lookup path should silently fall back to MC with an informative message.
  expect_message(
    val <- ExpectedClusteringInfoDistance(500L, nSim = 5L),
    "outside the shipped lookup range"
  )
  expect_true(is.finite(val) && val > 0 && val < 1)
  expect_message(
    ExpectedPhylogeneticInfoDistance(500L, nSim = 5L),
    "outside the shipped lookup range"
  )
})

test_that("ACID/APhID handle mixed-n tree lists", {
  trees <- list(
    BalancedTree(8L),
    PectinateTree(12L),
    BalancedTree(16L)
  )
  m_cid <- AdjustedClusteringInfoDistance(trees)
  m_pid <- AdjustedPhylogeneticInfoDistance(trees)
  # Returned as `dist` objects of length choose(3, 2) = 3.
  expect_s3_class(m_cid, "dist")
  expect_s3_class(m_pid, "dist")
  expect_length(m_cid, 3L)
  expect_length(m_pid, 3L)
  expect_true(all(is.finite(m_cid)))
  expect_true(all(is.finite(m_pid)))
  # Distance convention: 0 = identity, ~1 = random.  Worse-than-random
  # pairs can exceed 1, so we only sanity-check the lower bound here.
  expect_true(all(m_cid >= 0))
  expect_true(all(m_pid >= 0))
})

test_that("ExpectedRandomDistance rejects n < 4", {
  expect_error(ExpectedClusteringInfoDistance(3L), "must be an integer")
  expect_error(ExpectedPhylogeneticInfoDistance(3L), "must be an integer")
  expect_error(ExpectedClusteringInfoDistance(NA_integer_),
               "must be an integer")
})

test_that("`method = \"enumerate\"` is capped at n <= 7", {
  expect_error(ExpectedClusteringInfoDistance(8L, method = "enumerate"),
               "capped at n <= 7")
  expect_error(ExpectedPhylogeneticInfoDistance(8L, method = "enumerate"),
               "capped at n <= 7")
})

test_that("Random tree pairs at n = 20 average close to 1 under adjustment", {
  # A back-of-envelope sanity check: random pairs should score ~1 on
  # average, since E[D | random pair, n = 20] is exactly what we divide by.
  set.seed(2025)
  n <- 20L
  nReps <- 200L
  acid <- vapply(seq_len(nReps), function(i) {
    AdjustedClusteringInfoDistance(RandomTree(n, root = TRUE),
                                   RandomTree(n, root = TRUE))
  }, numeric(1))
  aphid <- vapply(seq_len(nReps), function(i) {
    AdjustedPhylogeneticInfoDistance(RandomTree(n, root = TRUE),
                                     RandomTree(n, root = TRUE))
  }, numeric(1))
  # Mean within ~3 SE of one (cid_sd / E / sqrt(nReps)).
  ref <- randomTreeDistances[randomTreeDistances[["n"]] == n, , drop = FALSE]
  se_cid <- (ref[["cid_sd"]] / ref[["cid_mean"]]) / sqrt(nReps)
  se_pid <- (ref[["pid_sd"]] / ref[["pid_mean"]]) / sqrt(nReps)
  expect_lt(abs(mean(acid)  - 1), 3 * se_cid)
  expect_lt(abs(mean(aphid) - 1), 3 * se_pid)
})

test_that("ACID/APhID handle non-overlapping tip sets", {
  # When two trees share no tips, CID/PID drop everything: the underlying
  # distance is well-defined but uninformative.  We just confirm the call
  # does not error and returns a finite value or NA.
  t1 <- BalancedTree(paste0("a", 1:8))
  t2 <- BalancedTree(paste0("b", 1:8))
  res_cid <- suppressWarnings(AdjustedClusteringInfoDistance(t1, t2))
  res_pid <- suppressWarnings(AdjustedPhylogeneticInfoDistance(t1, t2))
  expect_true(is.numeric(res_cid))
  expect_true(is.numeric(res_pid))
})

# -------------------------------------------------------------------------
# Shape-conditioned null tests (F4)
# -------------------------------------------------------------------------

test_that("shape null differs from size null at n=8 (pectinate vs balanced)", {
  t1 <- PectinateTree(8)
  t2 <- BalancedTree(8)
  acid_size  <- AdjustedClusteringInfoDistance(t1, t2, null = "size")
  set.seed(1)
  acid_shape <- AdjustedClusteringInfoDistance(t1, t2, null = "shape",
                                               method = "shape-mc",
                                               nSim = 1000L)
  # Shape null conditions on the specific shapes, so the expected distance
  # is different from the size null.  The two scores should not be equal.
  expect_false(isTRUE(all.equal(acid_size, acid_shape, tolerance = 1e-3)))

  aphid_size  <- AdjustedPhylogeneticInfoDistance(t1, t2, null = "size")
  set.seed(1)
  aphid_shape <- AdjustedPhylogeneticInfoDistance(t1, t2, null = "shape",
                                                  method = "shape-mc",
                                                  nSim = 1000L)
  expect_false(isTRUE(all.equal(aphid_size, aphid_shape, tolerance = 1e-3)))
})

test_that("shape-mc at n=6 converges to lookup within 4 SE", {
  # Shipped lookup currently covers n=4..7 (see data-raw/build_shape_lookup.R).
  # Use n=6 to exercise both shape-mc and shape-lookup paths.
  skip_if_not(
    nzchar(system.file("extdata", "shapeExpectedDistances.rds",
                       package = "TreeDist")),
    "shapeExpectedDistances.rds not built yet"
  )
  t1 <- PectinateTree(6)
  t2 <- BalancedTree(6)
  nSim <- 1e4L
  set.seed(8)
  mc_cid  <- ExpectedClusteringInfoDistance(null = "shape", method = "shape-mc",
                                            tree1 = t1, tree2 = t2, nSim = nSim)
  lk_cid  <- ExpectedClusteringInfoDistance(null = "shape", method = "shape-lookup",
                                            tree1 = t1, tree2 = t2)
  se_cid  <- attr(mc_cid, "sd") / sqrt(nSim)
  expect_equal(as.numeric(mc_cid), as.numeric(lk_cid),
               tolerance = 4 * se_cid)

  set.seed(8)
  mc_pid  <- ExpectedPhylogeneticInfoDistance(null = "shape", method = "shape-mc",
                                              tree1 = t1, tree2 = t2, nSim = nSim)
  lk_pid  <- ExpectedPhylogeneticInfoDistance(null = "shape", method = "shape-lookup",
                                              tree1 = t1, tree2 = t2)
  se_pid  <- attr(mc_pid, "sd") / sqrt(nSim)
  expect_equal(as.numeric(mc_pid), as.numeric(lk_pid),
               tolerance = 3 * se_pid)
})

test_that("shape null with shape-mc returns finite numeric at n=20 (beyond lookup)", {
  t1 <- PectinateTree(20)
  t2 <- BalancedTree(20)
  set.seed(20)
  # n=20 is beyond shipped lookup; should fall back to shape-mc gracefully
  val_cid <- suppressMessages(
    AdjustedClusteringInfoDistance(t1, t2, null = "shape",
                                   method = "shape-lookup", nSim = 200L)
  )
  expect_true(is.finite(val_cid))

  set.seed(20)
  val_cid_mc <- AdjustedClusteringInfoDistance(t1, t2, null = "shape",
                                               method = "shape-mc", nSim = 200L)
  expect_true(is.finite(val_cid_mc))
})

test_that("ACID(t, t, null = 'shape') == 0 (identity under shape null)", {
  for (n in c(6L, 8L, 12L)) {
    t1 <- PectinateTree(n)
    set.seed(n)
    val <- AdjustedClusteringInfoDistance(t1, t1, null = "shape",
                                          method = "shape-mc", nSim = 500L)
    expect_equal(val, 0,
                 label = paste0("ACID(pectinate(", n, "), same, shape)"))

    t2 <- BalancedTree(n)
    val2 <- AdjustedClusteringInfoDistance(t2, t2, null = "shape",
                                           method = "shape-mc", nSim = 500L)
    expect_equal(val2, 0,
                 label = paste0("ACID(balanced(", n, "), same, shape)"))
  }
})

test_that("shape null errors when trees are missing or have different n", {
  t1 <- PectinateTree(8)
  t2 <- BalancedTree(10)
  expect_error(
    ExpectedClusteringInfoDistance(null = "shape"),
    "must be supplied"
  )
  expect_error(
    ExpectedClusteringInfoDistance(null = "shape", tree1 = t1, tree2 = t2),
    "same number of tips"
  )
})
