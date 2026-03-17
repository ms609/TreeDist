## Tests for mutual_clustering_score() branches in src/pairwise_distances.cpp.
##
## These branches are only reachable via cpp_mutual_clustering_all_pairs(),
## invoked automatically by MutualClusteringInfo() when all trees share the
## same tip set and no R-level cluster is active (the fast path in
## .SplitDistanceAllPairs()).  The fast path returns a dist; MutualClusteringInfo
## then converts it to a full matrix and fills the diagonal with ClusteringEntropy
## values.  The off-diagonal entries equal the raw MCI scores from the batch C++
## function.

test_that("cpp_mutual_clustering_all_pairs: orthogonal splits score 0", {
  # Two 8-tip trees, each with exactly one internal split.
  # The splits cross orthogonally: every quadrant contains exactly 2 tips,
  # so a_and_b == A_and_b == a_and_B == A_and_B == 2.
  # This triggers the rounding-error guard (score = max_score), and
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

test_that("cpp_mutual_clustering_all_pairs: unequal split counts", {
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
  #       loop fills phantom rows ai=1,2 with max_score
  #       LAP solves the full 3×3 matrix
  #
  #   (col=1 → a=tr_3[3 splits], row=2 → b=tr_1b[1 split]):
  #     a_has_more = TRUE; most_splits = 3, b.n_splits = 1.
  #     For every ai=0..2: padRowAfterCol(ai, 1, max_score)

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

test_that("cpp_rf_info_all_pairs: multi-bin complement (> 64 tips)", {
  # Exercises b_complement[i][bin] = ~b.state[i][bin] at line 236, which
  # only runs when n_bins >= 2 (trees with > 64 tips on 64-bit platforms).
  skip_on_cran()
  trees <- ape::as.phylo(0:2, tipLabels = paste0("t", seq_len(100)))
  r <- InfoRobinsonFoulds(trees)
  expect_true(all(r >= 0))
  expect_equal(as.matrix(r)[2, 1],
               InfoRobinsonFoulds(trees[[1]], trees[[2]]),
               tolerance = 1e-10)
})

test_that("cpp_jaccard_all_pairs: allow_conflict = FALSE", {
  # Random trees have conflicting splits, so with allowConflict = FALSE
  # those pairs score max_score.
  trees <- ape::as.phylo(0:2, tipLabels = paste0("t", seq_len(20)))
  r <- JaccardRobinsonFoulds(trees, allowConflict = FALSE)
  expect_true(all(r >= 0))
  expect_equal(as.matrix(r)[2, 1],
               JaccardRobinsonFoulds(trees[[1]], trees[[2]],
                                     allowConflict = FALSE),
               tolerance = 1e-10)
})

test_that("cpp_jaccard_all_pairs: k = Inf and k != 1", {
  trees <- ape::as.phylo(0:2, tipLabels = paste0("t", seq_len(12)))
  # k = Inf triggers the isinf branch
  r_inf <- JaccardRobinsonFoulds(trees, k = Inf)
  expect_true(all(r_inf >= 0))
  expect_equal(as.matrix(r_inf)[2, 1],
               JaccardRobinsonFoulds(trees[[1]], trees[[2]], k = Inf),
               tolerance = 1e-10)
  # k = 2 triggers the pow() branch
  r_k2 <- JaccardRobinsonFoulds(trees, k = 2)
  expect_true(all(r_k2 >= 0))
  expect_equal(as.matrix(r_k2)[2, 1],
               JaccardRobinsonFoulds(trees[[1]], trees[[2]], k = 2),
               tolerance = 1e-10)
})

test_that("cpp_msd_all_pairs: unequal split counts, no exact matches", {
  # Three 6-tip trees with different numbers of non-trivial splits:
  #
  #   tr_3  — 3 splits: fully resolved
  #   tr_1a — 1 split:  {t1,t2,t3} | {t4,t5,t6}
  #   tr_1b — 1 split:  {t1,t4}    | {t2,t3,t5,t6}
  #
  # tr_1a and tr_1b share no splits (including complements), so their
  # comparison with tr_3 exercises:
  #   - line 361: padRowAfterCol when a.n_splits > b.n_splits
  #     (pair col=0 [tr_3, 3 splits] > row=1 [tr_1a, 1 split])
  #   - lines 393–396: the else branch (exact_n == 0) with unequal counts
  #     (pair col=0 [tr_3] vs row=1 [tr_1a], no shared splits)
  tr_3  <- ape::read.tree(text = "(((t1,t4),(t2,t5)),(t3,t6));")
  tr_1a <- ape::read.tree(text = "((t1,t2,t3),(t4,t5,t6));")
  tr_1b <- ape::read.tree(text = "((t1,t4),(t2,t3,t5,t6));")
  trees <- structure(list(tr_3, tr_1a, tr_1b), class = "multiPhylo")

  r <- MatchingSplitDistance(trees)
  m <- as.matrix(r)

  # Cross-validate batch vs single-pair path
  expect_equal(m[2, 1], MatchingSplitDistance(tr_3,  tr_1a), tolerance = 0)
  expect_equal(m[3, 1], MatchingSplitDistance(tr_3,  tr_1b), tolerance = 0)
  expect_equal(m[3, 2], MatchingSplitDistance(tr_1a, tr_1b), tolerance = 0)
})

test_that("Large trees: batch and single-pair paths agree (lg2 table bounds)", {
  # Exercise the lg2 / lg2_double_factorial / lg2_unrooted tables at indices

  # well beyond small-tree tests.  Catches out-of-bounds table access and
  # ensures the log2 decomposition in add_ic_element remains correct for
  # trees with many splits.
  skip_on_cran()
  set.seed(4728)
  trees <- ape::as.phylo(0:4, tipLabels = paste0("t", seq_len(200)))

  batch <- MutualClusteringInfo(trees)
  n <- length(trees)
  pairs <- which(upper.tri(batch), arr.ind = TRUE)
  for (k in seq_len(nrow(pairs))) {
    i <- pairs[k, 1]
    j <- pairs[k, 2]
    expect_equal(
      batch[i, j],
      MutualClusteringInfo(trees[[i]], trees[[j]]),
      tolerance = 1e-8,
      label = paste0("pair (", i, ",", j, ")")
    )
  }
})

test_that("Batch path handles asymmetric split counts with exact matches", {
  # Trees with different numbers of splits (polytomies) that share some

  # splits exactly.  After exact-match detection removes the shared splits,
  # the reduced LAP matrix is asymmetric (a_unmatched_n != b_unmatched_n),
  # exercising the padAfterRow guard in msd_score and jaccard_score.
  t1 <- ape::read.tree(text = "((t1,t2,t3),(t4,(t5,(t6,(t7,t8)))));")
  t2 <- ape::read.tree(text = "(((t1,t4),(t2,t3)),((t5,t6),(t7,t8)));")
  t3 <- ape::read.tree(text = "((t1,(t2,(t3,t4))),(t5,t6,t7,t8));")
  t4 <- ape::read.tree(text = "(((t1,t2),(t3,t4)),((t5,t7),(t6,t8)));")
  trees <- structure(list(t1, t2, t3, t4), class = "multiPhylo")

  # MSD batch vs per-pair
  msd_b <- as.matrix(MatchingSplitDistance(trees))
  for (i in seq_len(3)) for (j in (i + 1):4) {
    expect_equal(msd_b[i, j],
                 MatchingSplitDistance(trees[[i]], trees[[j]]),
                 tolerance = 1e-10,
                 label = paste0("MSD(", i, ",", j, ")"))
  }

  # Jaccard batch vs per-pair
  jrf_b <- as.matrix(JaccardRobinsonFoulds(trees))
  for (i in seq_len(3)) for (j in (i + 1):4) {
    expect_equal(jrf_b[i, j],
                 JaccardRobinsonFoulds(trees[[i]], trees[[j]]),
                 tolerance = 1e-10,
                 label = paste0("JRF(", i, ",", j, ")"))
  }
})

test_that("InfoRobinsonFoulds(similarity = TRUE) uses batch path", {
  # similarity = TRUE bypasses .FastDistPath() and falls through to
  # CalculateTreeDistance → .SplitDistanceAllPairs / .SplitDistanceManyMany,
  # exercising the IRF dispatch branches in tree_distance_utilities.R.
  trees <- ape::as.phylo(0:9, tipLabels = paste0("t", seq_len(20)))

  # All-pairs
  irf_sim <- InfoRobinsonFoulds(trees, similarity = TRUE)
  expect_true(inherits(irf_sim, "dist"))
  irf_mat <- as.matrix(irf_sim)
  expect_equal(irf_mat[2, 1],
               InfoRobinsonFoulds(trees[[1]], trees[[2]], similarity = TRUE),
               tolerance = 1e-10)

  # Cross-pairs
  tA <- trees[1:4]
  tB <- trees[5:10]
  irf_cross <- InfoRobinsonFoulds(tA, tB, similarity = TRUE)
  expect_equal(dim(irf_cross), c(4L, 6L))
  expect_equal(irf_cross[1, 1],
               InfoRobinsonFoulds(tA[[1]], tB[[1]], similarity = TRUE),
               tolerance = 1e-10)
})

test_that(".FastManyManyPath: cross-pairs distance exercises happy path", {
  # Exercises the happy path through .FastManyManyPath() in tree_distance.R:
  # cluster guard, TipLabels extraction, tip set matching, nTip check.
  tA <- ape::as.phylo(0:4, tipLabels = paste0("t", seq_len(20)))
  tB <- ape::as.phylo(5:9, tipLabels = paste0("t", seq_len(20)))

  cid_cross <- ClusteringInfoDistance(tA, tB)
  expect_equal(dim(cid_cross), c(5L, 5L))
  expect_equal(cid_cross[1, 1],
               ClusteringInfoDistance(tA[[1]], tB[[1]]),
               tolerance = 1e-10)
})

test_that(".FastManyManyPath: guards return NULL for edge cases", {
  tips <- paste0("t", seq_len(8))
  tA <- ape::as.phylo(0:2, tipLabels = tips)
  tB <- ape::as.phylo(3:5, tipLabels = tips)

  # Mismatched tip sets → falls back to slow path (returns non-NULL result)
  tB_diff <- ape::as.phylo(3:5, tipLabels = paste0("s", seq_len(8)))
  expect_true(!is.null(ClusteringInfoDistance(tA, tB_diff)))

  # Small trees (nTip < 4) → falls back to slow path
  tA3 <- ape::as.phylo(0:2, tipLabels = paste0("t", 1:3))
  tB3 <- ape::as.phylo(0:2, tipLabels = paste0("t", 1:3))
  expect_true(!is.null(ClusteringInfoDistance(tA3, tB3)))
})

test_that("SPI batch does not greedily match identical splits", {
  # Regression test for a correctness bug in shared_phylo_score().
  #
  # For the SPI (SharedPhylogeneticInfo) metric, matching identical splits
  # to each other is NOT necessarily globally optimal: spi_overlap(A, B)
  # where B *contains* A can exceed spi_overlap(A, A), so the full LAP
  # can find a better total score by pairing some identical splits with
  # non-identical (containing) partners.
  #
  # A greedy exact-match shortcut (match identical splits first, then
  # solve a reduced LAP) understates SPI, causing DifferentPhylogeneticInfo
  # (= Info1 + Info2 - 2 * SPI) to be overstated (trees appear more
  # different than they are).
  #
  # Counterexample: as.phylo(1, 10) and as.phylo(3, 10) share 5 of 7
  # splits.  The optimal LAP matches s1's 5|5 split to s2's 7|3 split
  # (which contains it) rather than to its identical copy in s2, because
  # spi_overlap({2,6,7,8,9,10}, {2,4,6,7,8,9,10}) > spi_overlap of the
  # identical pair.  The full-LAP SPI (32.012) exceeds the greedy
  # exact-match credit (31.771) even before the reduced LAP adds to it.
  # The same issue affects MatchingSplitInfoSplits (MSI).
  t1 <- ape::as.phylo(1, 10)
  t2 <- ape::as.phylo(3, 10)

  # Verify: optimal matching pairs s1[2] with s2[1] (non-identical),
  # even though s1[2] is identical to s2[2].
  s1 <- TreeTools::as.Splits(t1)
  s2 <- TreeTools::as.Splits(t2)
  result <- SharedPhylogeneticInfoSplits(s1, s2, reportMatching = TRUE)
  matching <- attr(result, "matching")

  # s1[2] is identical to s2[2] (same 5|5 split)
  expect_true(
    identical(as.logical(s1[[2]]), as.logical(s2[[2]])),
    label = "s1[2] and s2[2] are the same bipartition"
  )
  # But the optimal LAP matches s1[2] to s2[1] (a containing 7|3 split)
  expect_equal(matching[2], 1L,
               label = "optimal LAP matches s1[2] to s2[1], not its identical copy")

  # Batch path (all-pairs) must agree with the per-pair path
  trees <- structure(list(t1, t2), class = "multiPhylo")
  batch_spi <- SharedPhylogeneticInfo(trees)
  pair_spi <- SharedPhylogeneticInfo(t1, t2)
  expect_equal(as.matrix(batch_spi)[2, 1], pair_spi, tolerance = 1e-10)

  # The full-LAP SPI must exceed the sum of self-info for the 5 shared
  # splits (what the greedy approach would credit as exact-match score).
  # This proves the greedy is wrong — it understates SPI.
  n <- 10L
  shared_self_info <- sum(vapply(2:6, function(a) {
    k <- sum(as.logical(s1[[a]]))
    log2(TreeTools::NUnrooted(n)) -
      log2(TreeTools::NRooted(k)) - log2(TreeTools::NRooted(n - k))
  }, double(1)))
  expect_gt(pair_spi, shared_self_info,
            label = "full-LAP SPI exceeds greedy exact-match credit")

  # Verify the same holds for MSI (same vulnerability: msi_score formerly
  # used greedy exact-match detection, which is incorrect for MSI).
  batch_msi <- MatchingSplitInfoDistance(trees)
  pair_msi <- MatchingSplitInfoDistance(t1, t2)
  expect_equal(as.matrix(batch_msi)[2, 1], pair_msi, tolerance = 1e-10)
})
