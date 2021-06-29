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
                             p = c(ab = 0.8, d..x = 0.8, fgx = 0.6, 
                                   de = 0.8, cz = 0.6)),
               cons_phylo_info(trees))
  
  # Expected:
  # - tablei = 0: 
  #              1 2 3 4
  #   g..x : 2   x . . .
  #   f..x : 3   x x . .  ! 3|6, 60%
  #   d..e : 4   x x . x  ! 2|7, 80%
  #   d..x : 4   x x x .  ! 4|5, 80%
  #   z..x : 2   x . . .
  #   c..x : 4   x x x .  ! 7|2, 80%
  # 
  # - tablei = 1:
  #   Pass
  #   
  # - tablei = 2:
  #              3 4
  #   c..z : 3   x x  ! 2|7, 60%
  #   f..x : 1   . .  
  # 
  # - tablei = 3: Terminate.
  #
})
