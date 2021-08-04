library('TreeTools')

test_that("MAST fails gracefully", {
  expect_error(cpp_mast(BalancedTree(7)$edge, BalancedTree(8)$edge, 7)) # Different sizes
  expect_error(MASTSize(BalancedTree(8), UnrootTree(BalancedTree(8))))
  expect_error(MASTSize(BalancedTree(10000), PectinateTree(10000))) # Too large
  
})

test_that('MAST works', {
  tree1 <- BalancedTree(8L)
  tree2 <- PectinateTree(8L)
  expect_equal(8L, MASTSize(tree1, tree1, rooted = TRUE))
  expect_equal(8L, MASTSize(tree1, tree1, rooted = FALSE))
  expect_equal(4L, MASTSize(tree1, tree2, rooted = TRUE))
  expect_equal(6L, MASTSize(tree1, tree2, rooted = FALSE))
  
  expect_equal(MASTSize(BalancedTree(7), as.phylo(0:3, 7)),
               MASTSize(as.phylo(0:3, 7), BalancedTree(7)))
  
  expect_equal(MASTSize(list(BalancedTree(7), PectinateTree(7)), as.phylo(0:3, 7))[1, ],
               MASTSize(as.phylo(0:3, 7), BalancedTree(7)))
  
  expect_equal(MASTInfo(BalancedTree(7), as.phylo(0:3, 7)),
               MASTInfo(as.phylo(0:3, 7), BalancedTree(7)))
  
  expect_equal(MASTInfo(list(BalancedTree(7), PectinateTree(7)), as.phylo(0:3, 7))[1, ],
               MASTInfo(as.phylo(0:3, 7), BalancedTree(7)))
  
})

test_that("MAST supports funnily-ordered edges", {
  tree1 <- BalancedTree(8L)
  tree2 <- PectinateTree(8L)
  tree1$edge <- tree1$edge[c(1:7 * 2, (7:1 * 2) - 1), ]
  tree2$edge <- tree2$edge[c(1:7 * 2, (7:1 * 2) - 1), ]
  expect_equal(6L, MASTSize(Postorder(tree1), Postorder(tree2), rooted = FALSE))
  expect_equal(6L, MASTSize(tree1, tree2, rooted = FALSE))
})

test_that("MAST size calculated correctly on small trees", {
  library('TreeTools')
  #expect_equal(4L, MASTSize(as.phylo(0, 5), as.phylo(1, 5)))
  t1 <- structure(list(edge = structure(c(6L, 6L, 7L, 7L, 8L, 8L, 9L, 
                                          9L, 1L, 7L, 2L, 8L, 3L, 9L, 4L, 5L),
                                        .Dim = c(8L, 2L)), Nnode = 4L, 
                       tip.label = c("t3", "t5", "t4", "t1", "t2")),
                  class = "phylo", order = "cladewise")
  t2 <- structure(list(edge = structure(c(6L, 9L, 9L, 7L, 7L, 8L, 8L, 
                                          6L, 9L, 2L, 7L, 3L, 8L, 4L, 5L, 1L), 
                                        .Dim = c(8L, 2L)), 
                       tip.label = c("t3", "t4", "t1", "t2", "t5"), Nnode = 4L), 
                  class = "phylo", order = "cladewise")
  t1 <- RenumberTips(t1, paste0('t', 1:5))
  
  t2 <- RenumberTips(t2, t1)
  t2 <- Preorder(t2)

  ME <- function (e, node) {
    expect_equal(e, 
                 .MASTSizeEdges(Postorder(t1$edge),
                                RootOnNode(t2, node = node, TRUE)$edge,
                                nTip = 5L))
  }
  
  ME(2L, 1)
  ME(3L, 2)
  ME(4L, 3)
  ME(2L, 4)
  ME(2L, 5)
  ME(4L, 6)
  ME(4L, 7)
  ME(3L, 8)
  ME(3L, 9)
  
  expect_equal(4L, MASTSize(t1, t2, rooted = FALSE))
})
