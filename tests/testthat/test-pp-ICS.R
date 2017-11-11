context("Information content of steps")

test_that("Information content of steps calculated correctly", {
  expect_equal(c(1, 2), as.double(ICS(2, 2, 10000) * NUnrooted(4)))
  expect_equal(c(3, 12), signif(as.double(ICS(2, 3, 10000, warn=FALSE) * NUnrooted(5), 5)))
  
  expect_true(max(abs(cumsum(as.double(ICS(3,3, 6000, warn=FALSE) * NUnrooted(6))) - 
    c(NUnrootedMult(c(3,3)), NUnrootedMult(c(3,3)) + WithOneExtraStep(c(3,3)),
    NUnrooted(6)))) < 1)
    
  expect_true(max(abs(
    cumsum(as.double(ICS(3, 12, 60000, warn=FALSE))) - c(
      NUnrootedMult(c(3,12)) / NUnrooted(3 + 12),
      (NUnrootedMult(c(3,12)) + WithOneExtraStep(c(3,12))) / NUnrooted(3 + 12), 1)
  )) < 1e-03)
})