test_that("TransferDist() returns 0 for identical trees", {
  tr <- TreeTools::BalancedTree(8)
  expect_equal(0, TransferDist(tr, tr))
  expect_equal(0, TransferDist(tr, tr, scale = FALSE))
})


test_that("TransferDist() is symmetric", {
  bal8 <- TreeTools::BalancedTree(8)
  pec8 <- TreeTools::PectinateTree(8)
  expect_equal(TransferDist(bal8, pec8), TransferDist(pec8, bal8))
  expect_equal(TransferDist(bal8, pec8, scale = FALSE),
               TransferDist(pec8, bal8, scale = FALSE))
})


test_that("TransferDist() scaled is bounded by RF", {
  # The scaled transfer dissimilarity is bounded above by the RF distance
  bal8 <- TreeTools::BalancedTree(8)
  pec8 <- TreeTools::PectinateTree(8)
  td_scaled <- TransferDist(bal8, pec8)
  rf <- RobinsonFoulds(bal8, pec8)
  expect_lte(td_scaled, rf)
})


test_that("TransferDist() normalization", {
  bal8 <- TreeTools::BalancedTree(8)
  pec8 <- TreeTools::PectinateTree(8)

  # Scaled normalization denominator = 2 * (n - 3) = 10
  td_raw <- TransferDist(bal8, pec8)
  td_norm <- TransferDist(bal8, pec8, normalize = TRUE)
  expect_equal(td_norm, td_raw / (2 * (8 - 3)))
  expect_true(td_norm >= 0)
  expect_true(td_norm <= 1)

  # Unscaled normalization
  td_raw_u <- TransferDist(bal8, pec8, scale = FALSE)
  td_norm_u <- TransferDist(bal8, pec8, scale = FALSE, normalize = TRUE)
  # Even n: denom = (n^2 - 2n + 4) / 4 = (64 - 16 + 4) / 4 = 13
  expect_equal(td_norm_u, td_raw_u / 13)
})


test_that("TransferDist() all-pairs", {
  trees <- as.phylo(0:5, 8)

  # All-pairs
  d <- TransferDist(trees)
  expect_s3_class(d, "dist")
  expect_equal(attr(d, "Size"), 6L)

  # Check diagonal = 0 implicitly (dist objects don't store diagonal)
  # Check that batch matches pairwise
  dm <- as.matrix(d)
  for (i in 1:5) {
    for (j in (i + 1):6) {
      expected <- TransferDist(trees[[i]], trees[[j]])
      expect_equal(dm[i, j], expected, tolerance = 1e-12,
                   info = paste("pair", i, j))
    }
  }
})


test_that("TransferDist() cross-pairs", {
  trees1 <- as.phylo(0:2, 8)
  trees2 <- as.phylo(3:5, 8)

  mat <- TransferDist(trees1, trees2)
  expect_equal(dim(mat), c(3, 3))

  # Check matches pairwise
  for (i in 1:3) {
    for (j in 1:3) {
      expected <- TransferDist(trees1[[i]], trees2[[j]])
      expect_equal(mat[i, j], expected, tolerance = 1e-12,
                   info = paste("cross", i, j))
    }
  }
})


test_that("TransferDist() single tree vs list", {
  tr <- TreeTools::BalancedTree(8)
  trees <- as.phylo(0:3, 8)

  d1 <- TransferDist(tr, trees)
  d2 <- TransferDist(trees, tr)

  expect_length(d1, 4)
  expect_length(d2, 4)

  for (i in seq_along(trees)) {
    expected <- TransferDist(tr, trees[[i]])
    expect_equal(d1[i], expected, tolerance = 1e-12)
    expect_equal(d2[i], expected, tolerance = 1e-12)
  }
})


test_that("TransferDist() scaled vs unscaled relationship", {
  bal8 <- TreeTools::BalancedTree(8)
  pec8 <- TreeTools::PectinateTree(8)

  # Both should be non-negative
  td_s <- TransferDist(bal8, pec8, scale = TRUE)
  td_u <- TransferDist(bal8, pec8, scale = FALSE)
  expect_true(td_s >= 0)
  expect_true(td_u >= 0)

  # For identical trees, both = 0
  expect_equal(0, TransferDist(bal8, bal8, scale = TRUE))
  expect_equal(0, TransferDist(bal8, bal8, scale = FALSE))
})


test_that("TransferDist() reportMatching", {
  bal8 <- TreeTools::BalancedTree(8)
  pec8 <- TreeTools::PectinateTree(8)

  res <- TransferDist(bal8, pec8, reportMatching = TRUE)
  expect_equal(res[[1]], TransferDist(bal8, pec8))

  matching <- attr(res, "matching")
  expect_true(!is.null(matching))

  pairScores <- attr(res, "pairScores")
  expect_true(!is.null(pairScores))
  expect_equal(nrow(pairScores), TreeTools::NSplits(bal8))
  expect_equal(ncol(pairScores), TreeTools::NSplits(pec8))
})


test_that("TransferDist() handles star trees", {
  star <- TreeTools::StarTree(8)
  bal8 <- TreeTools::BalancedTree(8)

  # Star tree has no internal splits; dissimilarity = sum over binary tree only
  td <- TransferDist(bal8, star)
  expect_true(td >= 0)

  # Star vs star = 0
  expect_equal(0, TransferDist(star, star))
})


test_that("TransferDist() handles small trees", {
  # 4 tips: minimal non-trivial case
  tr1 <- TreeTools::BalancedTree(4)
  tr2 <- TreeTools::PectinateTree(4)

  td <- TransferDist(tr1, tr2)
  expect_true(is.finite(td))
  expect_true(td >= 0)
})


test_that("TransferDist() consistent with TransferConsensus internals", {
  skip_if_not_installed("TreeTools")

  # The transfer dissimilarity for a single pair should be consistent with
  # how TransferConsensus computes distances internally
  set.seed(6419)
  trees <- structure(lapply(1:10, function(i) ape::rtree(12)),
                     class = "multiPhylo")

  # All-pairs: check symmetry and non-negativity
  d <- TransferDist(trees)
  dm <- as.matrix(d)
  expect_true(all(dm >= -1e-12))
  expect_true(isSymmetric(dm, tol = 1e-12))
})


# =========================================================================
# Additional coverage: TransferDistance alias
# =========================================================================

test_that("TransferDistance() is an alias for TransferDist()", {
  bal8 <- TreeTools::BalancedTree(8)
  pec8 <- TreeTools::PectinateTree(8)
  expect_equal(TransferDistance(bal8, pec8), TransferDist(bal8, pec8))
  expect_equal(TransferDistance(bal8, pec8, scale = FALSE),
               TransferDist(bal8, pec8, scale = FALSE))
})


# =========================================================================
# cpp_transfer_dist() direct (both scaled + unscaled + matching)
# =========================================================================

test_that("cpp_transfer_dist() returns both scores and matchings", {
  bal8 <- TreeTools::BalancedTree(8)
  pec8 <- TreeTools::PectinateTree(8)
  tipLabels <- TreeTools::TipLabels(bal8)

  sp1 <- unclass(as.Splits(bal8, tipLabels))
  sp2 <- unclass(as.Splits(pec8, tipLabels))

  res <- TreeDist:::cpp_transfer_dist(sp1, sp2, 8L)

  expect_true(res$score_scaled >= 0)
  expect_true(res$score_unscaled >= 0)
  expect_length(res$matching_xy, nrow(sp1))
  expect_length(res$matching_yx, nrow(sp2))

  # Consistent with TransferDist
  expect_equal(res$score_scaled,
               TransferDist(bal8, pec8, scale = TRUE), tolerance = 1e-12)
  expect_equal(res$score_unscaled,
               TransferDist(bal8, pec8, scale = FALSE), tolerance = 1e-12)

  # Identical trees → distance 0, all matchings defined
  res_id <- TreeDist:::cpp_transfer_dist(sp1, sp1, 8L)
  expect_equal(res_id$score_scaled, 0)
  expect_equal(res_id$score_unscaled, 0)
})

test_that("cpp_transfer_dist() handles star trees", {
  star <- TreeTools::StarTree(8)
  bal8 <- TreeTools::BalancedTree(8)
  tipLabels <- TreeTools::TipLabels(bal8)

  sp_star <- unclass(as.Splits(star, tipLabels))
  sp_bal <- unclass(as.Splits(bal8, tipLabels))

  # Star vs binary: only one direction contributes
  res <- TreeDist:::cpp_transfer_dist(sp_bal, sp_star, 8L)
  expect_true(res$score_scaled >= 0)
  expect_true(res$score_unscaled >= 0)
  expect_length(res$matching_xy, nrow(sp_bal))
  # matching_yx has 0 elements for star tree
  expect_length(res$matching_yx, 0)

  # Star vs star: both scores = 0
  res_ss <- TreeDist:::cpp_transfer_dist(sp_star, sp_star, 8L)
  expect_equal(res_ss$score_scaled, 0)
  expect_equal(res_ss$score_unscaled, 0)
})

test_that("column mismatch caught by R-level guard", {
  sp1 <- matrix(as.raw(0), nrow = 1, ncol = 1)
  sp2 <- matrix(as.raw(0), nrow = 1, ncol = 2)
  expect_error(
    TreeDist:::.ValidateSplitArgs(sp1, sp2, 8L),
    "same number of tips"
  )
})


# =========================================================================
# cpp_transfer_dist_scored() edge cases
# =========================================================================

test_that("column mismatch caught by R-level guard (scored)", {
  sp1 <- matrix(as.raw(0), nrow = 1, ncol = 1)
  sp2 <- matrix(as.raw(0), nrow = 1, ncol = 2)
  expect_error(
    TreeDist:::.ValidateSplitArgs(sp1, sp2, 8L),
    "same number of tips"
  )
})

test_that("cpp_transfer_dist_scored() with asymmetric split counts", {
  bal8 <- TreeTools::BalancedTree(8)
  pec8 <- TreeTools::PectinateTree(8)
  tipLabels <- TreeTools::TipLabels(bal8)

  sp_bal <- unclass(as.Splits(bal8, tipLabels))
  sp_pec <- unclass(as.Splits(pec8, tipLabels))

  # BalancedTree(8) has fewer splits than PectinateTree(8)
  res <- TreeDist:::cpp_transfer_dist_scored(sp_bal, sp_pec, 8L, TRUE)
  expect_true(res$score >= 0)
  # matching length = max(nSplits1, nSplits2)
  expect_equal(length(res$matching),
               max(nrow(sp_bal), nrow(sp_pec)))
})


# =========================================================================
# Normalization edge cases
# =========================================================================

test_that(".TransferNormDenom() for odd and even nTip", {
  # Even nTip = 8: (64 - 16 + 4) / 4 = 13
  expect_equal(TreeDist:::.TransferNormDenom(8, FALSE), 13)
  # Odd nTip = 9: (81 - 18 + 5) / 4 = 17
  expect_equal(TreeDist:::.TransferNormDenom(9, FALSE), 17)
  # Scaled: 2 * (n - 3)
  expect_equal(TreeDist:::.TransferNormDenom(8, TRUE), 10)
  expect_equal(TreeDist:::.TransferNormDenom(9, TRUE), 12)
})

test_that("TransferDist() normalize with odd nTip", {
  trees <- as.phylo(0:3, nTip = 9)
  # Scaled normalization
  d <- TransferDist(trees, normalize = TRUE)
  d_raw <- TransferDist(trees, normalize = FALSE)
  denom <- 2 * (9 - 3)
  expect_equal(as.numeric(d), as.numeric(d_raw) / denom, tolerance = 1e-12)

  # Unscaled normalization with odd nTip
  d_u <- TransferDist(trees, scale = FALSE, normalize = TRUE)
  d_u_raw <- TransferDist(trees, scale = FALSE, normalize = FALSE)
  denom_u <- (81 - 18 + 5) / 4
  expect_equal(as.numeric(d_u), as.numeric(d_u_raw) / denom_u, tolerance = 1e-12)
})

test_that("TransferDist() normalize with cross-pairs", {
  trees1 <- as.phylo(0:2, 8)
  trees2 <- as.phylo(3:5, 8)

  mat_raw <- TransferDist(trees1, trees2)
  mat_norm <- TransferDist(trees1, trees2, normalize = TRUE)
  denom <- 2 * (8 - 3)
  expect_equal(mat_norm, mat_raw / denom, tolerance = 1e-12)

  # Unscaled cross-pairs normalize
  mat_raw_u <- TransferDist(trees1, trees2, scale = FALSE)
  mat_norm_u <- TransferDist(trees1, trees2, scale = FALSE, normalize = TRUE)
  denom_u <- 13
  expect_equal(mat_norm_u, mat_raw_u / denom_u, tolerance = 1e-12)
})

test_that("TransferDist() normalize via generic fallback (reportMatching)", {
  bal8 <- TreeTools::BalancedTree(8)
  pec8 <- TreeTools::PectinateTree(8)

  # reportMatching forces generic path; normalize divides after
  res <- TransferDist(bal8, pec8, normalize = TRUE, reportMatching = TRUE)
  raw <- TransferDist(bal8, pec8, reportMatching = TRUE)
  denom <- 2 * (8 - 3)
  expect_equal(res[[1]], raw[[1]] / denom, tolerance = 1e-12)

  # Unscaled + normalize + reportMatching
  res_u <- TransferDist(bal8, pec8, scale = FALSE,
                         normalize = TRUE, reportMatching = TRUE)
  raw_u <- TransferDist(bal8, pec8, scale = FALSE, reportMatching = TRUE)
  expect_equal(res_u[[1]], raw_u[[1]] / 13, tolerance = 1e-12)
})


# =========================================================================
# Fast path edge cases
# =========================================================================

test_that(".TransferDistAllPairs() returns NULL for edge cases", {
  # Single phylo (not multiPhylo)
  tr <- TreeTools::BalancedTree(8)
  expect_null(TreeDist:::.TransferDistAllPairs(tr, TRUE))

  # < 2 trees
  expect_null(TreeDist:::.TransferDistAllPairs(
    structure(list(tr), class = "multiPhylo"), TRUE))

  # < 4 tips
  small <- as.phylo(0:2, nTip = 3)
  expect_null(TreeDist:::.TransferDistAllPairs(small, TRUE))

  # Mismatched tip labels
  t1 <- TreeTools::BalancedTree(paste0("a", 1:6))
  t2 <- TreeTools::BalancedTree(paste0("b", 1:6))
  mixed <- structure(list(t1, t2), class = "multiPhylo")
  expect_null(TreeDist:::.TransferDistAllPairs(mixed, TRUE))
})

test_that(".TransferDistCrossPairs() returns NULL for edge cases", {
  tr1 <- TreeTools::BalancedTree(8)
  tr2 <- TreeTools::PectinateTree(8)

  # Both single → delegate to generic
  expect_null(TreeDist:::.TransferDistCrossPairs(tr1, tr2, TRUE))

  # < 4 tips
  small1 <- as.phylo(0:1, nTip = 3)
  small2 <- as.phylo(2:3, nTip = 3)
  expect_null(TreeDist:::.TransferDistCrossPairs(small1, small2, TRUE))

  # Mismatched tip labels
  t1 <- structure(list(TreeTools::BalancedTree(paste0("a", 1:6))),
                  class = "multiPhylo")
  t2 <- structure(list(TreeTools::BalancedTree(paste0("b", 1:6))),
                  class = "multiPhylo")
  expect_null(TreeDist:::.TransferDistCrossPairs(t1, t2, TRUE))
})

test_that("TransferDist() all-pairs with star trees in mix", {
  star <- TreeTools::StarTree(8)
  bal <- TreeTools::BalancedTree(8)
  pec <- TreeTools::PectinateTree(8)
  trees <- structure(list(star, bal, pec), class = "multiPhylo")

  d <- TransferDist(trees)
  dm <- as.matrix(d)
  expect_equal(dm[1, 1], 0)  # star self
  expect_true(dm[1, 2] >= 0) # star vs binary
  expect_equal(dm[2, 3], TransferDist(bal, pec), tolerance = 1e-12)
})

test_that("TransferDist() all-pairs with exactly 2 trees", {
  trees <- as.phylo(0:1, nTip = 8)
  d <- TransferDist(trees)
  expect_s3_class(d, "dist")
  expect_equal(attr(d, "Size"), 2L)
  expect_equal(as.numeric(d),
               TransferDist(trees[[1]], trees[[2]]), tolerance = 1e-12)
})

test_that("TransferDist() all-pairs with exactly 4 tips", {
  trees <- as.phylo(0:3, nTip = 4)
  d <- TransferDist(trees)
  expect_s3_class(d, "dist")
  dm <- as.matrix(d)
  expect_true(all(dm >= -1e-12))
})


# =========================================================================
# TransferDistSplits edge cases
# =========================================================================

test_that("TransferDistSplits() reportMatching with various configurations", {
  bal8 <- TreeTools::BalancedTree(8)
  pec8 <- TreeTools::PectinateTree(8)
  tipLabels <- TreeTools::TipLabels(bal8)

  sp_bal <- as.Splits(bal8, tipLabels)
  sp_pec <- as.Splits(pec8, tipLabels)

  # Unscaled reportMatching
  res <- TransferDistSplits(sp_bal, sp_pec, scale = FALSE,
                             reportMatching = TRUE)
  expect_true(is.numeric(res))
  expect_true(!is.null(attr(res, "matching")))
  expect_true(!is.null(attr(res, "pairScores")))
  expect_true(!is.null(attr(res, "matchedScores")))

  # Scaled reportMatching (already tested above, but confirm matchedSplits)
  res_s <- TransferDistSplits(sp_bal, sp_pec, scale = TRUE,
                               reportMatching = TRUE)
  expect_true(!is.null(attr(res_s, "matchedSplits")))
})

test_that("TransferDistSplits() with asymmetric nSplits (nSplits1 < nSplits2)", {
  # Use majority-rule consensus (fewer splits) vs fully resolved tree
  trees <- as.phylo(0:9, 8)
  mr <- Consensus(trees, p = 0.5)
  bal8 <- TreeTools::BalancedTree(8)
  tipLabels <- TreeTools::TipLabels(bal8)

  sp_mr <- as.Splits(mr, tipLabels)
  sp_bal <- as.Splits(bal8, tipLabels)

  # Confirm asymmetry
  expect_lt(nrow(sp_mr), nrow(sp_bal))

  res <- TransferDistSplits(sp_mr, sp_bal, reportMatching = TRUE)
  # Matching should be truncated to nSplits1
  matching <- attr(res, "matching")
  expect_equal(length(matching), nrow(sp_mr))
  expect_true(!is.null(attr(res, "pairScores")))
})

test_that("TransferDist() for 3-tip trees uses fallback path", {
  # 3 tips: fast path returns NULL, falls through to generic
  tr1 <- TreeTools::BalancedTree(3)
  tr2 <- TreeTools::PectinateTree(3)
  td <- TransferDist(tr1, tr2)
  expect_true(is.finite(td))
  expect_true(td >= 0)
})

test_that("TransferDist() with named tree list preserves labels", {
  trees <- as.phylo(0:3, nTip = 8)
  names(trees) <- paste0("tree_", 1:4)

  d <- TransferDist(trees)
  expect_equal(attr(d, "Labels"), paste0("tree_", 1:4))
})

test_that("TransferDist() cross-pairs with single Splits input", {
  bal8 <- TreeTools::BalancedTree(8)
  tipLabels <- TreeTools::TipLabels(bal8)
  sp <- as.Splits(bal8, tipLabels)
  trees <- as.phylo(0:3, nTip = 8)

  # Single Splits object vs list of trees
  result <- TransferDist(sp, trees)
  expect_length(result, 4)
  expect_true(all(result >= 0))
})

test_that("TransferDist fast-path guards handle edge cases", {
  # Tiny trees (< 4 tips) fall through to generic path
  tiny <- structure(rep(list(ape::read.tree(text = "(a,b,c);")), 3),
                    class = "multiPhylo")
  d <- TransferDist(tiny)
  expect_s3_class(d, "dist")

  # Cross-pairs with mixed-label trees fall through to generic path
  t1 <- ape::read.tree(text = "((a,b),(c,d));")
  t2 <- ape::read.tree(text = "((a,b),(e,f));")
  expect_true(is.numeric(TransferDist(t1, t2)))
})
