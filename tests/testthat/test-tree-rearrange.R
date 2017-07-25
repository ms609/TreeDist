library(ape)

context("tree rearrangement")
test_that("NNI works", {
  trComb    <- read.tree(text = "(((((1,2),3),4),5),6);")
  set.seed(0)
  nniComb <- NNI(trComb)
  expect_equal(nniComb$tip.label, trComb$tip.label)
  expect_equal(nniComb$Nnode, trComb$Nnode)
  
  
  trCombLen <- read.tree(text = "(((((1:1,2:1):1,3:2):1,4:3):1,5:1):1,6:1);")

  expect_equal(NNI(trComb), read.tree(text = "(((((3,2),1),4),5),6);"))
  
})
context("tree rearrangement")
test_that("TBR works", {
  tree <- read.tree(text="(((a, (b, (c, d))), (e, f)), (g, h));")
  library(phangorn)
  expect_equal(0, phangorn::RF.dist(TBR(tree, 3, 1), tree2 = read.tree(text="((a, ((b, (c, d)), (e, f))), (g, h));"), check.labels = TRUE))
  expect_warning(expect_identical(TBR(tree, 3, 2), tree))
  expect_warning(expect_identical(TBR(tree, 3, 3), tree))
  expect_warning(expect_identical(TBR(tree, 3, 2), tree))

  
  trCombLen <- read.tree(text = "(((((1:1,2:1):1,3:2):1,4:3):1,5:1):1,6:1);")

  expect_equal(NNI(trComb), read.tree(text = "(((((3,2),1),4),5),6);"))
  
})
