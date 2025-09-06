library("TreeTools")
test_that("list encoding tree works", {
  expect_equal(phylo_to_nested(BalancedTree(6)),
               list(list(list("t1", "t2"), "t3"),
                    list(list("t4", "t5"), "t6")))
  expect_equal(phylo_to_nested(PectinateTree(6)),
               list("t1", list("t2", list("t3", list("t4", list("t5", "t6"))))))
})
