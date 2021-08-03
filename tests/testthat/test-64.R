test_that("64-bit splits handled ok", {
  SPI <- function (nTip) SharedPhylogeneticInfo(TreeTools::BalancedTree(nTip), 
                                                TreeTools::PectinateTree(nTip),
                                                normalize = pmin)
  expect_gt(SPI(64), SPI(63))
  expect_lt(SPI(64), SPI(65))
  
  expect_gt(SPI(128), SPI(127))
  expect_lt(SPI(128), SPI(129))
})
