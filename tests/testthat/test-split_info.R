context("Day 1985: Split information")

test_that("Split info calculated", {
  library("ape")
  trees <- rev(list( # rev() because trees are emplaced backwards
    read.tree(text = '(a, (b, (c, (z, ((d, e), (f, (g, x)))))));'),
    read.tree(text = '(a, (b, (c, (z, ((d, e), (f, (g, x)))))));'),
    read.tree(text = '(a, (b, ((c, z), ((d, e), ((f, x), g)))));'),
    read.tree(text = '(a, (b, ((c, z), ((d, (e, x)), (f, g)))));'),
    read.tree(text = '(a, ((b, x), ((c, z), ((d, e), (f, g)))));')
  ))
  
  trees[] <- lapply(trees, RenumberTips, c(letters[1:7], 'x', 'z'))
  expect_equal(SplitwiseInfo(consensus(trees, p = 0.5), 
                             p = c(ab = 1, d..x = 0.8, fgx = 0.6, 
                                   de = 0.8, cz = 0.6)),
               cons_phylo_info(trees))
})
