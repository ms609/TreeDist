context("tree_distance_nni.R")
library('TreeTools')
library('ape')
test_that("Simple NNI approximations", {
  nTip <- 6L
  tree1 <- BalancedTree(nTip)
  tree2 <- PectinateTree(nTip)
  edge <- tree1$edge
  edge1 <- do.call(cbind,
                   PostorderEdges(edge[, 1], edge[, 2], dim(edge)[1], Nnode.phylo(tree1)))
  edge <- tree2$edge
  edge2 <- do.call(cbind,
                   PostorderEdges(edge[, 1], edge[, 2], dim(edge)[1], Nnode.phylo(tree2)))

  allMatched <- list(lower = 0L, tight_upper = 0L, loose_upper = 0L)
  oneUnmatched <- list(lower = 1L, tight_upper = 1L, loose_upper = 7L)
  fiveUnmatched <- list(lower = 5L, tight_upper = 10L, loose_upper = 17L + 4L)
  expect_equal(oneUnmatched, cpp_nni_distance(edge1, edge2, NTip(tree1)))
  expect_equal(oneUnmatched, NNIDist(tree1, tree2))

  # Identical trees
  tree1 <- Postorder(ape::read.tree(text="(((a, b), (c, d)), ((e, f), (g, h)));"))
  tree2 <- Postorder(ape::read.tree(text="(((a, b), (d, c)), ((h, g), (f, e)));"))
  expect_equal(allMatched, NNIDist(tree1, tree1))
  expect_equal(allMatched, NNIDist(tree1, tree2))
  
  # Only root edge is different
  tree2 <- Postorder(ape::read.tree(text="(((a, b), (e, f)), ((c, d), (g, h)));"))
  expect_equal(oneUnmatched, NNIDist(tree1, tree2))
  
  # Two separate regions of difference one
  tree2 <- ape::read.tree(text="((((a, b), c), d), (e, (f, (g, h))));")
  expect_equal(lapply(oneUnmatched, `*`, 2), NNIDist(tree1, tree2))
  
  # One region of three unmatched edges
  tree2 <- ape::read.tree(text="(((a, e), (c, d)), ((b, f), (g, h)));")
  expect_equal(list(lower = 3L, tight_upper = 5L, loose_upper = 13L), 
               NNIDist(tree1, tree2))
  
  # One region of four unmatched edges
  tree2 <- ape::read.tree(text="(((a, e), (f, d)), ((b, c), (g, h)));")
  expect_equal(list(lower = 4L, tight_upper = 7L, loose_upper = 18L), 
               NNIDist(tree1, tree2))
  
  # One region of five unmatched edges
  tree2 <- ape::read.tree(text="(((a, e), (f, d)), ((b, g), (c, h)));")
  expect_equal(fiveUnmatched, NNIDist(tree1, tree2))
  
  # Trees with different tips at root
  tree1 <- PectinateTree(1:8)
  tree2 <- ape::read.tree(text = '(3, ((5, 6), (7, (1, (2, (4, 8))))));')
  expect_equal(fiveUnmatched, NNIDist(tree1, tree2))
  
  # Too different for tight upper bound
  set.seed(10000)
  expect_true(is.na(NNIDist(ape::rtree(100, br=NULL), ape::rtree(100, br=NULL))[[2]]))
  
})

test_that("Crashes are avoided", {
  expect_error(NNIDist(list(PectinateTree(7), PectinateTree(8))))
  expect_error(NNIDist(list(PectinateTree(1:8), PectinateTree(8))))
  expect_error(NNIDist(list(PectinateTree(1:8), 
                            PectinateTree(as.character(1:8)))))
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
    lapply(NNIDist(list1), function (x) unname(as.matrix(x)[1:4, 4:1])),
    NNIDist(list1, rev(list1))
  )
})
  