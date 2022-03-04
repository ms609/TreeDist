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
  
  tree1 <- ape::read.tree(
    text = "(t3, ((t21, t17), (t4, (((((t9, t22), t25), t23), t24), t20))));")
  tree2 <- ape::read.tree(
    text = "(t3, ((t17, t21), (t4, (((((t9, t25), t23), t22), t24), t20))));")
  expect_true(all.equal(Reduce(tree1, tree2)[[1]],
                        DropTip(tree1, c("t3", "t17", "t21"))))
  
  liftRoot1 <- ape::read.tree(text = "(t1,(((((((((((((((((((((((((((((((((((t2,t70),t49),(t29,t63)),t80),(((t5,((t22,((t51,t113),t110)),t45)),t106),(t116,t123))),t109),t64),t69),t38),t73),t79),t28),(((t11,t40),t102),t33)),t18),t89),(t90,t117)),t34),t52),t9),((((t8,t127),(t91,t118)),t58),(t46,t48))),(((t12,t93),t74),t94)),(((((((t3,t4),t86),((((((((((((((t14,((t27,t54),(t39,(t53,t108)))),t119),(t17,((t25,((t32,(t95,t129)),t41)),t67))),(t19,t114)),t122),t42),t31),(((t15,t43),t84),t85)),((t66,t112),t104)),t71),(t44,(t47,t82))),t121),t65),(t56,(t92,t98)))),t81),(t16,t50)),((((((t6,t120),t60),t37),t105),t78),t23)),t75)),(t24,t126)),((((t124,t13),t20),t87),t68)),t103),t99),t59),(t10,((((t26,t61),t88),(t62,t115)),t30))),((t55,(t76,t96)),t107)),((t21,(t36,t77)),t101)),t130),t111),t72),t100),(t57,t128)));")
  liftRoot2 <- ape::read.tree(text =  "(t1,((((((((((((((((t2,t70),t102),(t11,t40)),t33),t28),t18),t89),t10),((((t26,t61),t88),(t62,t115)),t30)),((((((((((((((t3,t4),t86),((((((((((((((((((t5,((t22,((t51,t113),t110)),t45)),t106),(t116,t123)),(((t29,t63),t49),t80)),t109),t64),t38),t69),t79),(t73,t85)),(t43,t84)),(((((((t14,((t27,t54),(t39,(t53,t108)))),t119),(t17,((t25,((t32,t95),t41)),t67))),(t19,t114)),t122),t42),t31)),((t66,t112),t104)),t71),(t44,(t47,t82))),t121),t65),(t56,(t92,t98)))),t81),(t16,t50)),((((((t6,t120),t60),t37),t105),t78),t23)),t75),((((((t8,t127),(t91,t118)),(t13,t58)),(t46,t48)),(t9,((t15,t52),(t34,(t90,t117))))),((((t12,t93),t36),t74),t94))),t129),(t24,t126)),(((t124,t20),t87),t68)),t103),((t57,t128),t99)),t59)),((t55,(t76,t96)),t107)),((t21,t77),t101)),t130),t111),t72),t100));")
  expect_false("t70" %in% TipLabels(Reduce(liftRoot1, liftRoot2)))
  
})
