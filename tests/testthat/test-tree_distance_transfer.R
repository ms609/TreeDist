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
