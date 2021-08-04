test_that("Parallelization works", {
  library("TreeTools")
  trees <- as.phylo(0:20, 20)
  
  suppressMessages({
  expect_equal(options('TreeDist-cluster'), StartParallel(2))
  cl <- getOption('TreeDist-cluster')
  
  parallel <- ClusteringInfoDistance(trees)
  SetParallel(NULL)
  expect_null(GetParallel())
  expect_equal(ClusteringInfoDistance(trees), parallel)

  expect_equal(options('TreeDist-cluster'), SetParallel(cl))
  expect_equal(cl, getOption('TreeDist-cluster'))
  expect_equal(cl, GetParallel())
  parallel <- CompareAll(trees, NNIDist)
  SetParallel(NULL)
  expect_equal(CompareAll(trees, NNIDist), parallel)
  expect_false(StopParallel())
  SetParallel(cl)
  expect_true(StopParallel())
    })
})

