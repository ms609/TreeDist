test_that("k-means++ fails", {
  indistinct <- dist(rep(1, 100))
  expect_error(KMeansPP(indistinct), "Not enough distinct data points")
  expect_error(
    KMeansPP(as.matrix(indistinct)),
    "Not enough distinct data points"
  )
})

test_that("k-means++ works", {
  set.seed(0)
  two <- c(1:5, 11:15)
  expect_equal(KMeansPP(two, k = 1)$cluster, rep(1, 10))
  expect_equal(KMeansPP(two)$cluster, rep(1:2, each = 5))
  expect_equal(KMeansPP(cbind(numeric(10), two))$cluster,
                        rep(2:1, each = 5))
  five <- cbind(c(rnorm(10, -8), rnorm(5, 0), rnorm(10, 8)),
                c(rnorm(5, -5), rnorm(15, 5), rnorm(5, 0)))
  dists <- dist(five)
  expect_equal(unname(KMeansPP(dists, k = 1)$cluster), rep(1, 25))
  cl <- KMeansPP(dists, k = 5)$cluster
  # plot(five, col = KMeansPP(dists, k = 5)$cluster, pch = rep(15:19, each = 5))
  # plot(five, col = kmeans(dists, cent = 5)$cluster, pch = rep(15:19, each = 5))
  expect_equal(range(cl), c(1, 5))
  expect_equal(length(unique(cl[1:5])), 1)
  expect_equal(length(unique(cl[6:10])), 1)
  expect_equal(length(unique(cl[11:15])), 1)
  expect_equal(length(unique(cl[16:20])), 1)
  expect_equal(length(unique(cl[21:25])), 1)
})
