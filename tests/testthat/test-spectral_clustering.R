d <- dist(c(1:20, 60:41))

test_that("Spectral clustering fails gracefully", {
  expect_error(SpectralEigens(d, nEig = 0), "nEig must be.*positive")
})

test_that("Spectral clustering works", {
  allEig <- SpectralEigens(d, nEig = Inf)
  expect_equal(dim(allEig), c(40, 40))
  expect_equal(abs(SpectralEigens(d, nEig = 2)), abs(allEig[, 40:39]),
               tolerance = sqrt(.Machine[["double.eps"]]))

  # SpectralClustering is deprecated but should produce valid eigenvectors
  # Note: Column ordering may vary across BLAS implementations when
  # eigenvalues are equal, so we test properties rather than exact equality
  expect_warning(deprecated <- SpectralClustering(d, nEig = Inf),
                 "'SpectralClustering' is deprecated.")
  expect_equal(dim(deprecated), dim(allEig))

  # Verify eigenvectors are orthonormal (key property that should hold)
  eigen_prod <- t(deprecated) %*% deprecated
  expect_equal(eigen_prod, diag(ncol(deprecated)),
               tolerance = sqrt(.Machine[["double.eps"]]))
})
