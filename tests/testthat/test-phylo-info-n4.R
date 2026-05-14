# Regression test for the n=4 SharedPhylogeneticInfo / PhylogeneticInfoDistance
# bug discovered 2026-05-14.
#
# Root cause (src/pairwise_distances.cpp:713):
#   max_possible = lg2_unrooted[n_tips] - best_overlap
# At n_tips == 4, the only non-trivial split shape is balanced (2|2), so
# best_overlap == lg2_unrooted[4] == log2(3) and max_possible == 0.  The
# subsequent score_over_possible = max_score / 0 = +Inf produces NaN costs
# in the LAP, and shared_phylo_score returns 0 for every n=4 pair.
#
# Observable symptoms:
#   - SharedPhylogeneticInfo(t, t) returns 0 instead of SplitwiseInfo(t)
#   - PhylogeneticInfoDistance(t, t) returns 2 * SplitwiseInfo(t) instead of 0
#   - normalized PID is 1 for every n=4 pair, breaking ACID/APhID at n=4
#
# These tests are expected to FAIL until shared_phylo_score is patched to
# handle the max_possible == 0 edge case.  Once fixed, they pin the correct
# behaviour.

test_that("PhylogeneticInfoDistance(t, t) == 0 at n=4 (regression for max_possible bug)", {
  t <- ape::read.tree(text = "((a, b), (c, d));")
  expect_equal(PhylogeneticInfoDistance(t, t), 0)
  expect_equal(PhylogeneticInfoDistance(t, t, normalize = TRUE), 0)
})

test_that("SharedPhylogeneticInfo(t, t) == SplitwiseInfo(t) at n=4", {
  t <- ape::read.tree(text = "((a, b), (c, d));")
  expect_equal(SharedPhylogeneticInfo(t, t), SplitwiseInfo(t))
})

test_that("PhylogeneticInfoDistance is symmetric at n=4 for all 3 unrooted shapes", {
  trees <- list(
    ape::read.tree(text = "((a, b), (c, d));"),
    ape::read.tree(text = "((a, c), (b, d));"),
    ape::read.tree(text = "((a, d), (b, c));")
  )
  for (i in seq_along(trees)) {
    for (j in seq_along(trees)) {
      d_ij <- PhylogeneticInfoDistance(trees[[i]], trees[[j]])
      d_ji <- PhylogeneticInfoDistance(trees[[j]], trees[[i]])
      expect_equal(d_ij, d_ji)
      if (i == j) expect_equal(d_ij, 0)
    }
  }
})

test_that("PhylogeneticInfoDistance(t, t) == 0 at n=5..8 (control: should already pass)", {
  for (n in 5:8) {
    t <- TreeTools::BalancedTree(n)
    expect_equal(PhylogeneticInfoDistance(t, t), 0,
                 info = sprintf("n=%d", n))
    expect_equal(SharedPhylogeneticInfo(t, t), SplitwiseInfo(t),
                 info = sprintf("n=%d", n))
  }
})

test_that("ClusteringInfoDistance(t, t) == 0 at n=4 (control: CID branch unaffected)", {
  t <- ape::read.tree(text = "((a, b), (c, d));")
  expect_equal(ClusteringInfoDistance(t, t), 0)
  expect_equal(ClusteringInfoDistance(t, t, normalize = TRUE), 0)
})
