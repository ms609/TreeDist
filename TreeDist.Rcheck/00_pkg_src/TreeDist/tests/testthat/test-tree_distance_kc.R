library("TreeTools", quietly = TRUE)

test_that("KC fails gracefully", {
  expect_error(KendallColijn(BalancedTree(8), BalancedTree(1:8)),
               "Leaves must bear identical labels.")
})

test_that("KC vector calculations", {
  
  bal7 <- ape::read.tree(text = "(((t1,t2),(t3,t4)),((t5,t6),t7));")
  expect_error(PathVector(bal7$edge), "must be of class\\b")
  bal7b <- ape::read.tree(text = "(((t5,t6),t7), ((t1,t2),(t3,t4)));")
  expect_equal(PathVector(UnrootTree(bal7)),
               PathVector(UnrootTree(RootTree(bal7, 1))))
  expect_equal(PathVector(UnrootTree(bal7)),
               PathVector(UnrootTree(RenumberTips(bal7b, bal7))))
  expect_equal(as.numeric(PathVector(RenumberTips(bal7b, bal7))),
               c(2, 4, 4, 6, 6, 5,
                    4, 4, 6, 6, 5,
                       2, 6, 6, 5,
                          6, 6, 5,
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
  expect_equal(PathDist(trees), KendallColijn(trees, Vector = PathVector),
               ignore_attr = TRUE)
})

test_that("KCDiameter() calculated", {
  Test <- function(nTip) {
    tips <- seq_len(nTip)
    expect_equal(KendallColijn(PectinateTree(tips), PectinateTree(rev(tips))),
                 KCDiameter(nTip))
  }
  Test(4)
  Test(40)
  tree1 <- ape::read.tree(text = "(a, (b, (c, (d, (e, (f, (g, h)))))));")
  expect_equal(KendallColijn(tree1), 0)
  tree2 <- ape::read.tree(text = "(a, ((b, c), (d, (e, (f, (g, h))))));")
  tree3 <- ape::read.tree(text = "(a, (b, (c, (d, (e, (f, g))))));")
  trees <- c(tree1, tree2, tree3)
  expect_equal(KCDiameter(trees),
               c(KCDiameter(tree1), KCDiameter(tree2), KCDiameter(tree3)))
  expect_equal(KCDiameter(list(tree1, tree2, tree3)), KCDiameter(trees))
})
