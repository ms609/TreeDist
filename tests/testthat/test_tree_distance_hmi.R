library("TreeTools")
test_that("as.HPart works", {
  expect_equal(as.HPart(BalancedTree(6)),
               structure(class = "HPart",
                 list(list(list("t1", "t2"), list("t3")),
                      list(list("t4", "t5"), list("t6")))))
  expect_equal(as.HPart(PectinateTree(6)),
               structure(class = "HPart", 
                         list(list("t1"),
                              list(list("t2"),
                                   list(list("t3"),
                                        list(list("t4"),
                                             list("t5", "t6")))))))
})
