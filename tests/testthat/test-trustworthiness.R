library("TreeTools", quietly = TRUE)

test_that("MappingQuality() fails gracefully", {
  expect_error(ProjectionQuality(dist(1:10), dist(1:5)),
               "have the same dimensions")
})

test_that("Trustworthiness range", {
  
  for (N in c(5, 16, 50)) {
    k <- floor(seq(1, (N - 1) / 2, length.out = 5))
    UkMax <- pmin(k, N - k)
    maxVal <- N * ((UkMax * (N - k - 1)) - .Triangle(UkMax - 1))
    expect_equal(maxVal, (N * k * (2 * N - 3 * k - 1)) / 2)
  }
  
  trees <- as.phylo(1:12, 10)
  dists <- ClusteringInfoDist(trees)
  mapped <- dist(cmdscale(dists, k = 6))
  
  r <- apply(as.matrix(dists), 2, .Bercow) - 1
  diag(r) <- NA
  rHat <- apply(as.matrix(mapped), 2, .Bercow) - 1
  diag(r) <- NA
  
  rPrime <- 12 - r
  
  .Trustworthiness <- function (x, y, k) {
    MappingQuality(x, y, k)[["Trustworthiness"]]
  }
  
  expect_warning(expect_equal(.Trustworthiness(r, r, k = 12), NaN),
                 "All points are nearest neighbours")
  
  expect_equal(.Trustworthiness(r, r, k = 1), 1)
  expect_equal(.Trustworthiness(r, r, k = 5), 1)
  expect_equal(.Trustworthiness(r, r, k = 10), 1)
  
  expect_lte(.Trustworthiness(r, rHat, k = 1), 1)
  expect_lte(.Trustworthiness(r, rHat, k = 5), 1)
  expect_lte(.Trustworthiness(r, rHat, k = 10), 1)
  
  expect_gt(.Trustworthiness(r, rHat, k = 1), 0)
  expect_gt(.Trustworthiness(r, rHat, k = 5), 0)
  expect_gt(.Trustworthiness(r, rHat, k = 10), 0)
  
  expect_equal(.TrustSum(r, rPrime, 12, k = 1), .MMax(12, 1))
  expect_equal(.TrustSum(r, rPrime, 12, k = 5), .MMax(12, 5))
  expect_equal(.TrustSum(r, rPrime, 12, k = 10), .MMax(12, 10))
  
})
