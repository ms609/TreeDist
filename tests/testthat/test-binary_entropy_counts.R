test_that("binary_entropy_counts() fails gracefully", {
  expect_equal(binary_entropy_counts(integer(0), 0), numeric(0))
  expect_equal(binary_entropy_counts(integer(3), 0), rep(NA_real_, 3))
  expect_equal(binary_entropy_counts(integer(3), 2), rep(0, 3))
  expect_equal(binary_entropy_counts(c(1, NA, 2), 2), c(1, NA, 0))
})
