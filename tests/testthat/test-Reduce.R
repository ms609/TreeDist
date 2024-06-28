test_that("ReduceTrees() handles invalid input", {
  expect_error(ReduceTrees("STRING", as.phylo(1, 4)), "must be a `phylo`")
  expect_error(ReduceTrees(as.phylo(1, 4), logical(2)), "must be a `phylo`")
  expect_error(ReduceTrees(as.phylo(1, 4), as.phylo(1, 5)), "same leaf labels")
  expect_error(ReduceTrees(BalancedTree(5), StarTree(5)), "binary")
})

test_that("ReduceTrees()", {
  bal9 <- TreeTools::BalancedTree(9)
  pec9 <- TreeTools::PectinateTree(9)
  expect_null(ReduceTrees(bal9, bal9))
  par(mai = rep(0.1, 4), mfrow = c(2, 2))
  plot(RootTree(bal9, 1)); nodelabels()
  plot(RootTree(pec9, 1)); nodelabels()
  confl <- ReduceTrees(bal9, pec9)
  expect_true(
    all.equal(confl[[1]],
              ape::read.tree(text = "(t1, ((t4, t5), ((t6, t7), t9)));")))
  expect_true(all.equal(confl[[2]], DropTip(pec9, c(2, 3, 8))))
  
  pec9b <- RenumberTips(TreeTools::PectinateTree(paste0('t', c(2:9, 1))),
                        pec9)
  pecred <- ReduceTrees(pec9b, pec9)
  expect_true(all.equal(pecred[[2]], PectinateTree(5)))
  expect_true(all.equal(pecred[[1]], PectinateTree(paste0('t', c(1, 5:2)))))
  
  long1 <- ape::read.tree(
    text = "(a, (b, (c, (d, (e, (f, ((g, (X, Y)), (h, (i, j)))))))));")
  long2 <- ape::read.tree(
    text = "(b, (c, (d, (e, (f, (g, ((h, (X, Y)), (i, (j, a)))))))));")
  longRed <- ReduceTrees(long1, long2)
  expect_true(all.equal(longRed[[1]],
                        DropTip(long1, c("b", "c", "X"))))
  expect_true(all.equal(longRed[[2]],
                        RootTree(DropTip(long2, c("b", "c", "X")), "a")))
  
  long1 <- ape::read.tree(
    text = "(r, (oo, (t, (a, (b, (c, (d, (e, (f, ((g, (X, Y)), (h, (i, j))))))))))));")
  long2 <- ape::read.tree(
    text = "(r, (oo, (t, (b, (c, (d, (e, (f, (g, ((h, (X, Y)), (i, (j, a))))))))))));")
  longRed <- ReduceTrees(long1, long2)
  expect_true(all.equal(longRed[[1]],
                        DropTip(long1, c("t", "oo", "b", "c", "X"))))
  expect_true(all.equal(longRed[[2]],
                        DropTip(long2, c("t", "oo", "b", "c", "X"))))
  
  tree1 <- ape::read.tree(
    text = "(t3, ((t21, t17), (t4, (((((t9, t22), t25), t23), t24), t20))));")
  tree2 <- ape::read.tree(
    text = "(t3, ((t17, t21), (t4, (((((t9, t25), t23), t22), t24), t20))));")
  expect_true(all.equal(ReduceTrees(tree1, tree2)[[1]],
                        DropTip(tree1, c("t4", "t17", "t21", "t20", "t24"))))

  # Lift root AND reduce chain to new root base
  keptRoot1 <- ape::read.tree(text = "((a,(b,((((c,(d,(((e,(((f,((g,(h,((((((i,j),k),l),m),n),o))),p)),q),r)),s),t))),u),v),x))),root);")
  keptRoot2 <- ape::read.tree(text = "((a,(b,(((c,(d,(((e,(((f,((g,(h,(((((((x,i),j),k),l),m),n),o))),p)),q),r)),s),t))),u),v))),root);")
  tree1 <- RenumberTips(keptRoot1, c("root", letters))
  tree2 <- keptRoot2
  expect_equal(c("root", "x", "a", "b") %in% TipLabels(ReduceTrees(tree1, tree2)),
               c(TRUE, TRUE, FALSE, FALSE))
})
