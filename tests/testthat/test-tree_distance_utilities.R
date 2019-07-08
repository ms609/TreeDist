context('tree_distance_utilities.R')
test_that('Tree normalization works', {
  expect_equal(0.5, NormalizeInfo(5, 3, 5, how=TRUE, 
                                  InfoInTree=function(x, y) x + y, y = 1L))
  expect_equal(0:4 / c(1, 2, 3, 3, 3), NormalizeInfo(0:4, 1:5, 3, InfoInTree = I,
                                       Combine = min))
  expect_equal(3/4, NormalizeInfo(unnormalized = 3, 1, 1, how = 4L))
  expect_equal(matrix(0:24 / 10, 5, 5), 
               NormalizeInfo(0:24, rep(5, 5), rep(5, 5), InfoInTree = I))
  expect_equal(matrix(1:20 / 10, 5, 4) * c(10,10,1,1,1, 10,10,1,1,1, rep(1,10)), 
               NormalizeInfo(1:20, c(1, 1, 10, 10, 10), c(1, 1, 10, 10), 
                             InfoInTree = I, Combine = pmax))
})
