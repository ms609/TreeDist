test_that("PathDist()", {
  library("TreeTools", quietly = TRUE, warn.conflicts = FALSE)
  expect_equal(c(5.66, 6, 6, 6.32, 6.32, 5.74),
               PathDist(as.phylo(0:5, 6), BalancedTree(paste0("t", 6:1))),
               tolerance = 2)
  expect_equal(c(5.66, 6, 6, 6.32, 6.32, 5.74),
               PathDist(BalancedTree(6), as.phylo(0:5, 6)),
               tolerance = 2)
  
  expect_equal(PathDist(BalancedTree(6), PectinateTree(6)),
               PathDist(list(BalancedTree(6), PectinateTree(6)))[1],
               ignore_attr = TRUE)
  
  trees <- as.phylo(1:8, 29)
  expect_equal(unname(as.matrix(PathDist(trees))), PathDist(trees, trees))
  expect_equal(PathDist(trees), CompareAll(trees, PathDist),
               ignore_attr = TRUE)
})

test_that("PathDist() equivalent to path.dist()", {
  skip_if_not_installed("phangorn")
  postTrees <- Postorder(as.phylo(0:5, 182))
  expect_equal(PathDist(postTrees), phangorn::path.dist(postTrees))
  ub <- microbenchmark::microbenchmark
  ub(PathDist(postTrees), phangorn::path.dist(postTrees))
  pv <- profvis::profvis
  pv(ub(PathDist(postTrees), phangorn::path.dist(postTrees)))
})
  