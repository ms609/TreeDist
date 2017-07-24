library(ape)

context("tree search")
test_that("NNI works", {
  trComb    <- read.tree(text = "(((((1,2),3),4),5),6);")
  trCombNNI <- 10
  trCombLen <- read.tree(text = "(((((1:1,2:1):1,3:2):1,4:3):1,5:1):1,6:1);")

  set.seed(0)
  expect_equal(NNI(trComb, dataset), )
  
})
