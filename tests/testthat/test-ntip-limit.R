test_that(".SL_MAX_TIPS is populated", {
  expect_true(is.integer(TreeDist:::.SL_MAX_TIPS))
  expect_true(TreeDist:::.SL_MAX_TIPS >= 2048L)
})

test_that("Trees exceeding SL_MAX_TIPS are rejected", {
  limit <- TreeDist:::.SL_MAX_TIPS
  # Build a mock Splits raw matrix with limit + 1 tips.
  # The matrix needs: nrow = limit - 2 splits, ncol = ceil((limit + 1) / 8)
  bad_nTip <- limit + 1L
  n_splits <- bad_nTip - 3L
  n_cols <- ceiling(bad_nTip / 8)
  mock <- matrix(as.raw(0), nrow = n_splits, ncol = n_cols)
  class(mock) <- c("Splits", class(mock))
  attr(mock, "nTip") <- bad_nTip
  attr(mock, "tip.label") <- paste0("t", seq_len(bad_nTip))

  expect_error(
    TreeDist:::GeneralizedRF(mock, mock, bad_nTip,
                             TreeDist:::cpp_robinson_foulds_distance,
                             maximize = FALSE, reportMatching = FALSE),
    "exceed"
  )
})
