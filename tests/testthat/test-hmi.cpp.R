library("TreeTools")

test_that("HMI calculated correctly", {
  bal6 <- BalancedTree(6)
  pec6 <- PectinateTree(6)
  expect_equal(d_n_nested(phylo_to_nested(bal6),
                          phylo_to_nested(pec6)),
               HierachicalMutual(bal6, pec6))
})
