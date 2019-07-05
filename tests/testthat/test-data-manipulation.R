context("data_manipulation.R")
test_that("minimum step counts are correctly calculated", {
  expect_equal(1, MinimumSteps(1:3))
  expect_equal(1, MinimumSteps(c(1:3, 5)))
  expect_equal(0, MinimumSteps(c(6, 7, 14)))
  expect_equal(1, MinimumSteps(0:3)) # 0 representing the inapplicable token

  
  # 04, 14, 24, 34, 05, 16, 27, 38, 9A
  # In this case, chosing the most common state (4) means that we have to choose 567&8 too
  # 012&3 is a better solution
  # We also have to choose one of 9 or A, but it doesn't matter which.
  expect_equal(4, MinimumSteps(c(
    2^0 + 2^4,
    2^1 + 2^4,
    2^2 + 2^4,
    2^3 + 2^4,
    2^0 + 2^5,
    2^1 + 2^6,
    2^2 + 2^7,
    2^3 + 2^8,
    2^9 + 2^10
  )))
  
  
})
