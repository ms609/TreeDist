context("Partitions.R")

test_that("UniqueSplits works", {
  suppressWarnings(RNGversion("3.5.0")) # Until we can require R3.6.0
  set.seed(1)
  splits6 <- Tree2Splits(ape::rtree(6, br=NULL))
  expect_equal(c('8'=FALSE, '10'=FALSE, '11'=FALSE), UniqueSplits(splits6)['t4', ])
  expect_equal(!splits6, UniqueSplits(cbind(!splits6, splits6), TRUE))
  
})

test_that("Large splits don't cause memory issues", {
  splits5000 <- cbind(c(rep(TRUE, 2), rep(FALSE, 4998)),
                      c(rep(FALSE, 2), rep(TRUE, 4998)),
                      c(rep(TRUE, 2500), rep(FALSE, 2500)),
                      c(rep(TRUE, 2500), rep(FALSE, 2500)),
                      c(rep(TRUE, 2500), rep(FALSE, 2500)),
                      c(rep(FALSE, 2500), rep(TRUE, 2500))
  )
  expect_equal(c(5000, 2), dim(UniqueSplits(splits5000)))
})
