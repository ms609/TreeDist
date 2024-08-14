d <- dist(c(1:20, 60:41))

test_that("Spectral clustering fails gracefully", {
  expect_error(SpectralEigens(d, nEig = 0), "nEig must be.*positive")
})

test_that("Spectral clustering works", {
  allEig <- SpectralEigens(d, nEig = Inf)
  expect_equal(dim(allEig), c(40, 40))
  expect_equal(SpectralEigens(d, nEig = 2), allEig[, 40:39])
})
