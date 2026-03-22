test_that("Header mutual_clustering_score matches MutualClusteringInfo", {
  skip_if_not_installed("TreeDist")
  library(TreeTools)
  library(TreeDist)

  bal8 <- BalancedTree(8)
  pec8 <- PectinateTree(8)
  star8 <- StarTree(8)

  tips <- TipLabels(bal8)
  n_tip <- length(tips)
  splits_bal <- as.Splits(bal8, tips)
  splits_pec <- as.Splits(pec8, tips)
  splits_star <- as.Splits(star8, tips)

  # Score-only from the installable-header implementation
  impl_score <- TreeDist:::cpp_mci_impl_score

  impl_bal_pec <- impl_score(splits_bal, splits_pec, n_tip)
  impl_bal_bal <- impl_score(splits_bal, splits_bal, n_tip)
  impl_star <- impl_score(splits_bal, splits_star, n_tip)

  # Reference from MutualClusteringInfo (unnormalized score)
  ref_bal_pec <- MutualClusteringInfo(bal8, pec8)
  ref_bal_bal <- MutualClusteringInfo(bal8, bal8)

  expect_equal(impl_bal_pec, ref_bal_pec, tolerance = 1e-10)
  expect_equal(impl_bal_bal, ref_bal_bal, tolerance = 1e-10)
  expect_equal(impl_star, 0)
})

test_that("Header MCI covers exact-match early exit and partial LAP", {
  skip_if_not_installed("TreeDist")
  library(TreeTools)
  library(TreeDist)
  impl_score <- TreeDist:::cpp_mci_impl_score

  # Two identical trees → all splits match exactly (early exit path)
  bal20 <- BalancedTree(20)
  tips <- TipLabels(bal20)
  n_tip <- length(tips)
  splits20 <- as.Splits(bal20, tips)

  result <- impl_score(splits20, splits20, n_tip)
  expect_equal(result, MutualClusteringInfo(bal20, bal20), tolerance = 1e-10)

  # Trees that share some but not all splits → partial match + LAP
  pec20 <- PectinateTree(20)
  splits_pec20 <- as.Splits(pec20, tips)

  result2 <- impl_score(splits20, splits_pec20, n_tip)
  expect_equal(result2, MutualClusteringInfo(bal20, pec20), tolerance = 1e-10)
})
