test_that("KC vector calculations", {
  bal7 <- ape::read.tree(text = "(((t1,t2),(t3,t4)),((t5,t6),t7));")
  bal7b <- ape::read.tree(text = "(((t5,t6),t7), ((t1,t2),(t3,t4)));")
  expect_equal(PathVector(bal7), PathVector(RootTree(bal7, 1)))
  expect_equal(PathVector(bal7), PathVector(bal7b))
  expect_equal(as.numeric(PathVector(bal7b)),
               c(2, 4, 4, 5, 5, 4,
                    4, 4, 5, 5, 4,
                       2, 5, 5, 4,
                          5, 5, 4,
                             2, 3,
                                3))
  
  expect_equal(SplitVector(bal7), SplitVector(RootTree(bal7, 1)))
  expect_equal(SplitVector(bal7), SplitVector(bal7b))
  expect_equal(as.numeric(SplitVector(BalancedTree((9)))),
               c(2, 3, 5, 5, 7, 7, 7, 7,
                    3, 5, 5, 7, 7, 7, 7,
                       5, 5, 7, 7, 7, 7,
                          2, 6, 6, 6, 6,
                             6, 6, 6, 6,
                                2, 4, 4,
                                   4, 4,
                                      2))
})

test_that("KC distances with special vectors", {
  trees <- as.phylo(1:20, 12)
  expect_equivalent(PathDist(trees), KendallColijn(trees, Vector = PathVector))
})
