test_that("k-means++ works", {
  set.seed(0)
  expect_equal(KMeansPP(1:10, k = 1)$cluster, rep(1, 10))
  expect_equal(KMeansPP(1:10)$cluster, rep(1:2, each = 5))
  expect_equal(KMeansPP(cbind(1:10, rep(1:2, each = 5)))$cluster,
                        rep(2:1, each = 5))
})
