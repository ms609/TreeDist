test_that("Split info calculated", {
  library("ape")
  trees <- rev(list( # rev() because trees are emplaced backwards
    read.tree(text = '(a, (b, (c, (z, ((d, e), (f, (g, x)))))));'),
    read.tree(text = '(a, (b, (c, (z, ((d, e), (f, (g, x)))))));'),
    read.tree(text = '(a, (b, ((c, z), ((d, e), ((f, x), g)))));'),
    read.tree(text = '(a, (b, ((c, z), ((d, (e, x)), (f, g)))));'),
    read.tree(text = '(a, ((b, x), ((c, z), ((d, e), (f, g)))));')
  ))
  split_p <- c(ab = 0.8, d..x = 0.8, fgx = 0.6, de = 0.8, cz = 0.6)
  
  trees[] <- lapply(trees, RenumberTips, c(letters[1:7], 'x', 'z'))
  expect_equal(consensus_info(trees, TRUE, 1), 0)
  expect_equal(SplitwiseInfo(consensus(trees, p = 0.5), p = split_p),
  expect_warning(expect_equal(ConsensusInfo(trees, p = 0.4),
                              ConsensusInfo(trees, p = 0.5)))
  expect_equal(SplitwiseInfo(consensus(trees, p = 0.7),
                             p = split_p[split_p > 0.7]),
               consensus_info(trees, TRUE, 0.7))
  expect_error(consensus_info(trees, TRUE, 7))
  expect_error(consensus_info(trees, TRUE, -7))
  
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
  
  p_in <- c(7, 5, 3, 2, 2) / 9
  p_out <- 1 - p_in
  expect_equal(sum(apply(rbind(p_in, p_out), 2, Entropy) * split_p)
               * NTip(trees[[1]]),
               consensus_info(trees, FALSE, 0.5))
  
  # Even number of trees: cz, with p == 0.5, not in consensus.
  split_p <- c(ab = 0.8, d..x = 0.8, fgx = 0.6, de = 0.8, cz = 0.6)
  expect_equal(SplitwiseInfo(consensus(trees[-1], p = 0.5),
                             p = c(ab = 1, d..x = 1, fgx = 3/4, de = 3/4)),
               consensus_info(trees[-1], TRUE, p = 0.5))
  expect_equal(SplitwiseInfo(consensus(trees[-1], p = 0.8),
                             p = c(ab = 1, d..x = 1)),
               consensus_info(trees[-1], TRUE, p = 0.8))
  expect_warning(expect_equal(ConsensusInfo(trees[-1], p = 80),
                              ConsensusInfo(trees[-1], p = 1)))
})
