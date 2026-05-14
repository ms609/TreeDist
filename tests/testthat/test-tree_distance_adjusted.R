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

test_that("ACID/APhID return 1 for identical trees", {
  for (n in c(6L, 20L, 50L)) {
    trees <- list(
      balanced = BalancedTree(n),
      pectinate = PectinateTree(n)
    )
    set.seed(n)
    trees[["random"]] <- RandomTree(n, root = TRUE)
    for (tr in trees) {
      expect_equal(AdjustedClusteringInfoDistance(tr, tr), 1)
      expect_equal(AdjustedPhylogeneticInfoDistance(tr, tr), 1)
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
  # n = 300 is beyond the shipped table's max (n = 200); the lookup path
  # should silently fall back to MC with an informative message.
  expect_message(
    val <- ExpectedClusteringInfoDistance(300L, nSim = 5L),
    "outside the shipped lookup range"
  )
  expect_true(is.finite(val) && val > 0 && val < 1)
  expect_message(
    ExpectedPhylogeneticInfoDistance(300L, nSim = 5L),
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
  expect_true(all(m_cid <= 1))
  expect_true(all(m_pid <= 1))
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

test_that("Random tree pairs at n = 20 average close to 0 under adjustment", {
  # A back-of-envelope sanity check: random pairs should score ~0 on
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
  # Mean within ~3 SE of zero (cid_sd / E / sqrt(nReps)).
  ref <- randomTreeDistances[randomTreeDistances[["n"]] == n, , drop = FALSE]
  se_cid <- (ref[["cid_sd"]] / ref[["cid_mean"]]) / sqrt(nReps)
  se_pid <- (ref[["pid_sd"]] / ref[["pid_mean"]]) / sqrt(nReps)
  expect_lt(abs(mean(acid)), 3 * se_cid)
  expect_lt(abs(mean(aphid)), 3 * se_pid)
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
