library(ape)
library(phangorn)

context("tree search")
comb11 <- read.tree(text="(a, (b, (c, (d, (e, (f, (g, (h, (i, (j, k))))))))));")
data11 <- cbind(upper.tri(matrix(FALSE, 11, 11))[, 3:10], lower.tri(matrix(FALSE, 11, 11))[, 2:9])
rownames(data11) <- letters[1:11]
phy11 <- phyDat(data11, type='USER', levels=c(FALSE, TRUE))

test_that("Tree can be found", {
  set.seed(0)
  expect_equal(DoTreeSearch(RandomTree(phy11, 'a'), phy11, maxIter=2500, Rearrange = RootedTBR, verbosity=0), comb11)
  expect_equal(DoTreeSearch(RandomTree(phy11, 'a'), phy11, maxIter=2000, Rearrange = RootedSPR, verbosity=0), comb11)
  expect_equal(DoTreeSearch(RandomTree(phy11, 'a'), phy11, maxIter=2000, Rearrange = RootedNNI, verbosity=0), comb11)
  expect_equal(Ratchet(RandomTree(phy11, 'a'), phy11, searchIter=300, searchHits = 20, ratchHits=3, verbosity=-1), comb11)
  expect_equal(SectorialSearch(RandomTree(phy11, 'a'), phy11, searchIter=300, searchHits = 20, ratchHits=3, verbosity=-1), comb11)
})
