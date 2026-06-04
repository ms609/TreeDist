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

test_that("k-means++ matrix seeding matches materialized distance matrix", {
  # On-the-fly distance rows must select the same centres as the original
  # implementation, which materialized the full n * n distance matrix.  The
  # oracle below is the literal pre-refactor body of KMeansPP.matrix(); both run
  # under the same seed, so the comparison checks that the refactor preserved
  # the RNG draw sequence and chose identical centres (hence cluster and score).
  oracle <- function(x, k, nstart) {
    n <- dim(x)[[1]]
    ret <- list(tot.withinss = Inf)
    d <- as.matrix(dist(x))
    for (start in seq_len(nstart)) {
      centres <- integer(k)
      centres[1L] <- sample.int(n, 1L)
      min_d <- d[centres[1L], ]
      for (i in 2:k) {
        centres[i] <- sample.int(n, 1L, prob = min_d ^ 2)
        min_d <- pmin.int(min_d, d[centres[i], ])
      }
      proposal <- kmeans(x, centers = x[centres, ], iter.max = 100L)
      if (proposal[["tot.withinss"]] < ret[["tot.withinss"]]) {
        ret <- proposal
      }
    }
    ret
  }
  set.seed(1)
  x <- matrix(rnorm(500 * 5), ncol = 5)
  set.seed(42)
  expected <- oracle(x, k = 8, nstart = 5)
  set.seed(42)
  got <- KMeansPP(x, k = 8, nstart = 5, iter.max = 100L)
  expect_identical(got[["cluster"]], expected[["cluster"]])
  expect_identical(got[["tot.withinss"]], expected[["tot.withinss"]])
})
