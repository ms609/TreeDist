library(ape)
library(phangorn)

context("tree search")
comb11 <- read.tree(text="(a, (b, (c, (d, (e, (f, (g, (h, (i, (j, k))))))))));")
data11 <- cbind(upper.tri(matrix(FALSE, 11, 11))[, 3:10], lower.tri(matrix(FALSE, 11, 11))[, 2:9])
rownames(data11) <- letters[1:11]
phy11 <- phyDat(data11, type='USER', levels=c(FALSE, TRUE))

test_that("Tree can be found", {
  set.seed(0)
  expect_equal(Ratchet(RandomTree(phy11), phy11, searchIter=1000, ratchHits=6, outgroup='1'), comb11)
})
