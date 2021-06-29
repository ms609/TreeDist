context("Day 1985: Split information")

test_that("Split info calculated", {
  library("ape")
  trees <- list(
    read.tree(text = '(a, (b, (c, (d, e), (f, (g, x)))));'),
    read.tree(text = '(a, (b, (c, (d, e), (f, (g, x)))));'),
    read.tree(text = '(a, (b, (c, (d, e), ((f, x), g))));'),
    read.tree(text = '(a, (b, (c, (d, (e, x)), (f, g))));'),
    read.tree(text = '(a, ((b, x), (c, (d, e), (f, g))));')
  )
  expect_equal(SplitwiseInfo(consensus(trees, p = 0.5), p = c(1, 0.6, 0.8)),
               cons_phylo_info(trees))
})
