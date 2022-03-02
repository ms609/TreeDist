test_that("Reduce()", {
  tree1 <- TreeTools::BalancedTree(9)
  tree2 <- TreeTools::PectinateTree(9)
  expect_null(Reduce(tree1, tree1))
  par(mai = rep(0.1, 4), mfrow=c(2, 2))
  plot(RootTree(tree1, 1)); nodelabels()
  plot(RootTree(tree2, 1)); nodelabels()
  confl <- Reduce(tree1, tree2)
  expect_true(
    all.equal(confl[[1]],
              ape::read.tree(text = "(t1, ((t4, t5), ((t6, t7), t8)));")))
  expect_true(all.equal(confl[[2]], DropTip(tree2, c(2, 3, 9))))
})
