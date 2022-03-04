test_that("Reduce()", {
  bal9 <- TreeTools::BalancedTree(9)
  pec9 <- TreeTools::PectinateTree(9)
  expect_null(Reduce(bal9, bal9))
  par(mai = rep(0.1, 4), mfrow=c(2, 2))
  plot(RootTree(bal9, 1)); nodelabels()
  plot(RootTree(pec9, 1)); nodelabels()
  confl <- Reduce(bal9, pec9)
  expect_true(
    all.equal(confl[[1]],
              ape::read.tree(text = "(t1, ((t4, t5), ((t6, t7), t9)));")))
  expect_true(all.equal(confl[[2]], DropTip(pec9, c(2, 3, 8))))
  
  pec9b <- RenumberTips(TreeTools::PectinateTree(paste0('t', c(2:9, 1))),
                        pec9)
  plot(RootTree(pec9b, 1)); nodelabels()
  plot(RootTree(pec9, 1)); nodelabels()
  pec9b$edge[PostorderOrder(pec9b), ]
  pec9$edge[PostorderOrder(pec9), ]
  pecred <- Reduce(pec9b, pec9)
  expect_true(all.equal(pecred[[2]], PectinateTree(5)))
  expect_true(all.equal(pecred[[1]], PectinateTree(paste0('t', c(1, 5:2)))))
  
  long1 <- ape::read.tree(
    text = "(a, (b, (c, (d, (e, (f, ((g, (X, Y)), (h, (i, j)))))))));")
  long2 <- ape::read.tree(
    text = "(b, (c, (d, (e, (f, (g, ((h, (X, Y)), (i, (j, a)))))))));")
  longRed <- Reduce(long1, long2)
  expect_true(all.equal(longRed[[1]],
                        DropTip(long1, c("b", "c", "X"))))
  expect_true(all.equal(longRed[[2]],
                        RootTree(DropTip(long2, c("b", "c", "X")), "a")))
  
  long1 <- ape::read.tree(
    text = "(r, (oo, (t, (a, (b, (c, (d, (e, (f, ((g, (X, Y)), (h, (i, j))))))))))));")
  long2 <- ape::read.tree(
    text = "(r, (oo, (t, (b, (c, (d, (e, (f, (g, ((h, (X, Y)), (i, (j, a))))))))))));")
  longRed <- Reduce(long1, long2)
  expect_true(all.equal(longRed[[1]],
                        DropTip(long1, c("b", "oo", "t", "c", "X"))))
  expect_true(all.equal(longRed[[2]],
                        DropTip(long2, c("b", "oo", "t", "c", "X"))))
})
