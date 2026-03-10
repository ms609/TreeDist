## Tests for mutual_clustering_score() branches in src/pairwise_distances.cpp.
##
## These branches are only reachable via cpp_mutual_clustering_all_pairs(),
## invoked automatically by MutualClusteringInfo() when all trees share the
## same tip set and no R-level cluster is active (the fast path in
## .SplitDistanceAllPairs()).  The fast path returns a dist; MutualClusteringInfo
## then converts it to a full matrix and fills the diagonal with ClusteringEntropy
## values.  The off-diagonal entries equal the raw MCI scores from the batch C++
## function.

test_that("cpp_mutual_clustering_all_pairs: orthogonal splits score 0 (line 82)", {
  # Two 8-tip trees, each with exactly one internal split.
  # The splits cross orthogonally: every quadrant contains exactly 2 tips,
  # so a_and_b == A_and_b == a_and_B == A_and_B == 2.
  # This triggers the rounding-error guard (score = max_score) at line 82, and
  # the LAP assigns that sole pair, yielding MCI = 0.
  tr1 <- ape::read.tree(text = "((t1,t2,t3,t4),(t5,t6,t7,t8));")
  tr2 <- ape::read.tree(text = "((t1,t2,t5,t6),(t3,t4,t7,t8));")
  trees <- structure(list(tr1, tr2), class = "multiPhylo")

  r <- MutualClusteringInfo(trees)

  # r is a 2×2 matrix; the off-diagonal [2,1] is the pairwise MCI
  expect_equal(r[2, 1], 0, tolerance = 1e-10)
  # Agrees with the single-pair path (cpp_mutual_clustering, not the batch fn)
  expect_equal(r[2, 1], MutualClusteringInfo(tr1, tr2), tolerance = 1e-10)
})

test_that("cpp_mutual_clustering_all_pairs: unequal split counts (lines 94, 131-138)", {
  # Three 6-tip trees with different numbers of non-trivial splits:
  #
  #   tr_1a — 1 split:  {t1,t2,t3} | {t4,t5,t6}
  #   tr_3  — 3 splits: {t1,t4}, {t2,t5}, {t3,t6}
  #   tr_1b — 1 split:  {t1,t2}   | {t3,t4,t5,t6}
  #
  # The batch loop (col < row) processes three pairs using 0-indexed trees:
  #
  #   (col=0 → a=tr_1a[1 split], row=1 → b=tr_3[3 splits]):
  #     a_has_more = FALSE; most_splits = 3.
  #     No split pairs are exact matches, so exact_n = 0 → else branch:
  #       loop fills phantom rows ai=1,2 with max_score        (lines 131–134)
  #       LAP solves the full 3×3 matrix                       (lines 136–138)
  #
  #   (col=1 → a=tr_3[3 splits], row=2 → b=tr_1b[1 split]):
  #     a_has_more = TRUE; most_splits = 3, b.n_splits = 1.
  #     For every ai=0..2: padRowAfterCol(ai, 1, max_score)    (line 94)

  tr_1a <- ape::read.tree(text = "((t1,t2,t3),(t4,t5,t6));")
  tr_3  <- ape::read.tree(text = "(((t1,t4),(t2,t5)),(t3,t6));")
  tr_1b <- ape::read.tree(text = "((t1,t2),(t3,t4,t5,t6));")
  trees <- structure(list(tr_1a, tr_3, tr_1b), class = "multiPhylo")

  r <- MutualClusteringInfo(trees)

  # r is a 3×3 matrix; off-diagonal [i, j] with i > j is pair (j, i)
  expect_equal(dim(r), c(3L, 3L))

  # Cross-validate all three off-diagonal pairs against the single-pair path
  # (which uses cpp_mutual_clustering in tree_distances.cpp, not the batch fn)
  expect_equal(r[2, 1], MutualClusteringInfo(tr_1a, tr_3),  tolerance = 1e-9)
  expect_equal(r[3, 1], MutualClusteringInfo(tr_1a, tr_1b), tolerance = 1e-9)
  expect_equal(r[3, 2], MutualClusteringInfo(tr_3,  tr_1b), tolerance = 1e-9)
})
