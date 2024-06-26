test_that("PathDist()", {
  library("TreeTools", quietly = TRUE)
  t05 <- as.phylo(0:5, 6)
  bal6 <- BalancedTree(paste0("t", 6:1))
  vec6 <- PathVector(RenumberTips(bal6, t05))
  expect_equal(PathDist(t05, bal6),
               sqrt(colSums((vapply(t05, PathVector, vec6) - vec6) ^ 2)))
  expect_equal(PathDist(bal6, t05), PathDist(t05, bal6))
  
  expect_equal(PathDist(bal6, PectinateTree(6)),
               PathDist(list(bal6, PectinateTree(6)))[1],
               ignore_attr = TRUE)
  
  trees <- as.phylo(1:8, 29)
  expect_equal(unname(as.matrix(PathDist(trees))), PathDist(trees, trees))
  expect_equal(PathDist(trees), CompareAll(trees, PathDist),
               ignore_attr = TRUE)
})

test_that("PathDist() equivalent to path.dist()", {
  skip_if_not_installed("phangorn")
  
  t1 <- Postorder(as.phylo(0:5, 6))
  t2 <- Postorder(BalancedTree(6))
  t3 <- Postorder(PectinateTree(6))
  expect_equal(PathDist(UnrootTree(t1), UnrootTree(t2)),
               phangorn::path.dist(t1, t2))
  expect_equal(PathDist(UnrootTree(t2), UnrootTree(t3)),
               phangorn::path.dist(t3, t2))
  
  postTrees <- Postorder(as.phylo(0:5, 182))
  expect_equal(PathDist(postTrees), phangorn::path.dist(postTrees))
})
  