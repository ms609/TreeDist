library("TreeTools")
cl <- makeCluster(2)
SetParallel(cl)
trees <- as.phylo(0:20, 20)

test_that("TreeDist parallelizes", {
  SetParallel(cl)
  parallel <- ClusteringInfoDistance(trees)
  SetParallel(NULL)
  on.exit(SetParallel(cl))
  expect_equal(ClusteringInfoDistance(trees), parallel)
})

test_that("CompareAll() parallelizes", {
  SetParallel(cl)
  parallel <- CompareAll(trees, NNIDist)
  SetParallel(NULL)
  on.exit(SetParallel(cl))
  expect_equal(CompareAll(trees, NNIDist), parallel)
})

StopParallel()
