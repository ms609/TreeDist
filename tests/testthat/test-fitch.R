library(ape)
library(phangorn)

context("tree scoring")
comb11 <- read.tree(text="(a, (b, (c, (d, (e, (f, (g, (h, (i, (j, k))))))))));")
data11 <- cbind(upper.tri(matrix(FALSE, 11, 11))[, 3:10], lower.tri(matrix(FALSE, 11, 11))[, 2:9])
rownames(data11) <- letters[1:11]
phy11 <- phyDat(data11, type='USER', levels=c(FALSE, TRUE))

test_that("Tree is scored correctly", {
  expect_equal(FitchScore(comb11, phy11), 16)
  expect_equal(FitchScore(Pruningwise(comb11), phy11), 16)
  expect_equal(FitchScore(Cladewise(comb11), phy11), 16)
  
  expect_equal(FitchScore(comb11, data11, TipData=TipsAreRows), 16)
})
