context("tree_distance_nni.R")
library('TreeTools')

test_that("NNIDist() handles exceptions", {
  expect_error(NNIDist(list(PectinateTree(7), PectinateTree(8))))
  expect_error(NNIDist(list(PectinateTree(1:8), PectinateTree(8))))
  expect_error(NNIDist(list(PectinateTree(1:8), 
                            PectinateTree(as.character(1:8)))))
  expect_error(cpp_nni_distance( # Too many tips.
    PectinateTree(40000)$edge, # Will fail before not being postorder is problem
    BalancedTree(40000)$edge, 40000))
  
  expect_error(NNIDist(BalancedTree(5), RootOnNode(BalancedTree(5), 1)))
  
})

test_that("Simple NNI approximations", {
  nTip <- 6L
  tree1 <- BalancedTree(nTip)
  tree2 <- PectinateTree(nTip)
  edge1 <- Postorder(tree1$edge)
  edge2 <- Postorder(tree2$edge)
  
  allMatched <- c(lower = 0L, tight_upper = 0L, loose_upper = 0L)
  oneUnmatched <- c(lower = 1L, tight_upper = 1L, loose_upper = 7L)
  fiveUnmatched <- c(lower = 5L, tight_upper = 10L, loose_upper = 17L + 4L)
  
  Test <- function (expectation, tree) {
    expect_equal(expectation, NNIDist(tree1, tree))
    for (i in c(2, 3, 4, 6)) {
      tree1i <- RootOnNode(tree1, i)
      j <- 0
      for (t2 in unique(lapply(1:9, RootOnNode, tree = tree))) {
        expect_equal(expectation, NNIDist(tree1i, t2))
      }
    }
  }
  
  expect_equal(allMatched, NNIDist(BalancedTree(2), PectinateTree(2)))
  
  expect_equal(oneUnmatched, cpp_nni_distance(edge1, edge2, NTip(tree1)))
  Test(oneUnmatched, PectinateTree(nTip))

  # Identical trees
  tree1 <- Postorder(ape::read.tree(text="(((a, b), (c, d)), ((e, f), (g, h)));"))
  tree2 <- Postorder(ape::read.tree(text="(((a, b), (d, c)), ((h, g), (f, e)));"))
  Test(allMatched, tree1)
  Test(allMatched, tree2)
  
  # Tree names
  output <- NNIDist(list(bal = tree1, pec = tree2), 
                    as.phylo(0:2, tipLabels = letters[1:8]))
  expect_equal(rownames(output), c('bal', 'pec'))
  
  # Only root edge is different
  Test(oneUnmatched, 
       Postorder(ape::read.tree(text="(((a, b), (e, f)), ((c, d), (g, h)));")))
  
  # Two separate regions of difference one
  Test(oneUnmatched * 2, 
       ape::read.tree(text="((((a, b), c), d), (e, (f, (g, h))));"))
  
  # One region of three unmatched edges
  Test(c(lower = 3L, tight_upper = 5L, loose_upper = 13L), 
       ape::read.tree(text="(((a, e), (c, d)), ((b, f), (g, h)));"))
  
  # One region of four unmatched edges
  Test(c(lower = 4L, tight_upper = 7L, loose_upper = 18L),
       tree2 <- ape::read.tree(text="(((a, e), (f, d)), ((b, c), (g, h)));"))
  
  # One region of five unmatched edges
  Test(fiveUnmatched, 
       ape::read.tree(text="(((a, e), (f, d)), ((b, g), (c, h)));"))
  
  # Trees with different leaves at root
  tree1 <- PectinateTree(1:8)
  Test(fiveUnmatched, 
       ape::read.tree(text = '(3, ((5, 6), (7, (1, (2, (4, 8))))));'))
  
  # Too different for tight upper bound
  set.seed(10000)
  expect_true(is.na(NNIDist(rtree(100, br = NULL), rtree(100, br = NULL))[[2]]))
  
})

test_that("NNI with lists of trees", {
  tree1 <- BalancedTree(1:8)
  list1 <- list(tree1, PectinateTree(as.character(1:8)),
                PectinateTree(as.character(c(4:1, 5:8))), BalancedTree(c(1:3, 8:4)))
  
  multResult <- NNIDist(tree1, list1)
  expect_equal(NNIDist(tree1, list1[[1]]), multResult[, 1])
  expect_equal(NNIDist(tree1, list1[[2]]), multResult[, 2])
  expect_equal(NNIDist(tree1, list1[[3]]), multResult[, 3])
  expect_equal(NNIDist(tree1, list1[[4]]), multResult[, 4])
  
  expect_equal(NNIDist(tree1, list1), NNIDist(list1, tree1))
  
  # CompareAll
  expect_equal(CompareAll(list1, NNIDist), NNIDist(list1))
  
  expect_equivalent(
    vapply(NNIDist(list1), function (x) unname(as.matrix(x)[1:4, 4:1]),
           matrix(0,4,4)),
    NNIDist(list1, rev(list1))
  )
})
  