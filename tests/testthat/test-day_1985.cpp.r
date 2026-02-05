test_that("Day 1985 does not overflow", {
  bigTree <- PectinateTree(2^14 + 1)
  expect_no_error(TreeTools::as.ClusterTable(bigTree))
  expect_equal(RobinsonFoulds(list(bigTree, bigTree))[[1]], 0)
})

test_that("Day 1985 examples", {
  library("TreeTools", quietly = TRUE)
  
  PrepareTree <- function(text) {
    tmp <- ape::read.tree(text = text)
    RenumberTips(tmp, as.character(seq_along(tmp$tip.label)))
  }
  
  TestRF <- function(t1, t2) {
    expect_equal(RobinsonFoulds(t1, t2), 
                 NSplits(t1) + NSplits(t2) - (2 * COMCLUST(list(t1, t2))))
  }
  
  t1 <- PrepareTree("(1, (2, (3, (4, (5, 6)))));")
  t2 <- PrepareTree("(1, (2, ((4, 3), (6, 5))));")
  t3 <- PrepareTree("(1, (((2, 3), 5), (4, 6)));")
  tStar <- PrepareTree("(1, 2, 3, 4, 5, 6);")
  c1 <- as.ClusterTable(t1)
  c2 <- as.ClusterTable(t2)
  c3 <- as.ClusterTable(t3)
  cStar <- as.ClusterTable(tStar)
  t1p <- ape::read.tree(text = write.tree(t1))
  t2p <- ape::read.tree(text = write.tree(t2))
  tStarp <- ape::read.tree(text = write.tree(tStar))
  
  expect_equal(2L, robinson_foulds_all_pairs(list(c1, c2)))
  
  TestRF(t1, t2)
  TestRF(t1, t3)
  TestRF(t2, t3)
  TestRF(t1, tStar)
  TestRF(tStar, tStar)
  expect_equal(c(RobinsonFoulds(t1, c(t2, t3, tStar)), RobinsonFoulds(t2, t3),
                 NSplits(c(t2, t3))),
               as.integer(RobinsonFoulds(c(t1, t2, t3, tStar))))
  
  # Day's figure 2
  t1 <- PrepareTree("((10, 7), (6, (8, 11)), (12, (4, (2, 1))), 14, (5, 9, 13), 3);")
  as.ClusterTable(t1)
  t2 <- PrepareTree("(((2, 4, 5, 7, 9, 10, 12, 13), 1, 14), (6, (8, 11)), 3);")
  as.ClusterTable(t2)
  
  expect_equal(2L, COMCLUST(list(t1, t2)))
  expect_equal(7L, as.integer(RobinsonFoulds(list(t1, t2))))

})
  
test_that("RobinsonFoulds() with realistic trees", {
  # Test exposes segfault that has arisen whilst modernizing C++
  tr1 <- structure(list(
    edge = structure(c(8L, 8L, 9L, 10L, 10L, 9L, 11L, 11L, 8L, 12L, 12L,
                       1L, 9L, 10L, 2L, 3L, 11L, 4L, 5L, 12L, 6L, 7L),
                     dim = c(11L, 2L)),
    Nnode = 5L, tip.label = c("t1", "t2", "t3", "t4", "t5", "t6", "t7")),
    class = "phylo", order = "preorder")
  tr2 <- structure(list(
    edge = structure(c(8L, 9L, 10L, 10L, 9L, 11L, 11L, 8L, 12L, 12L, 8L,
                       9L, 10L, 1L, 2L, 11L, 3L, 4L, 12L, 5L, 6L, 7L),
                     dim = c(11L, 2L)),
    Nnode = 5L, tip.label = c("t1", "t2", "t3", "t4", "t5", "t6", "t7")),
    class = "phylo", order = "preorder")
  t4 <- list(a = tr1, b = tr2, c = tr1, d = tr2)
  r4 <- TreeTools::RootTree(t4, 1)
  
  expect_equal(as.numeric(RobinsonFoulds(r4)),
               c(8, 0, 8, 8, 0, 8))
})
