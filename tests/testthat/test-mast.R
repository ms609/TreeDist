library("TreeTools")

test_that("MAST fails gracefully", {
  # R-level guard catches edge count mismatch
  expect_error(.MASTSizeEdges(BalancedTree(7)$edge, BalancedTree(8)$edge, 7),
               "same number of edges")
  expect_error(MASTSize(BalancedTree(8), UnrootTree(BalancedTree(8))))
  expect_error(MASTSize(BalancedTree(10000), PectinateTree(10000))) # Too large
  
})

test_that("MAST works", {
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
  library("TreeTools")
  #expect_equal(4L, MASTSize(as.phylo(0, 5), as.phylo(1, 5)))
  t1 <- structure(list(edge = structure(c(6L, 6L, 7L, 7L, 8L, 8L, 9L,
                                          9L, 1L, 7L, 2L, 8L, 3L, 9L, 4L, 5L),
                                        dim = c(8L, 2L)), Nnode = 4L,
                       tip.label = c("t3", "t5", "t4", "t1", "t2")),
                  class = "phylo", order = "cladewise")
  t2 <- structure(list(edge = structure(c(6L, 9L, 9L, 7L, 7L, 8L, 8L,
                                          6L, 9L, 2L, 7L, 3L, 8L, 4L, 5L, 1L),
                                        dim = c(8L, 2L)), 
                       tip.label = c("t3", "t4", "t1", "t2", "t5"), Nnode = 4L), 
                  class = "phylo", order = "cladewise")
  t1 <- RenumberTips(t1, paste0("t", 1:5))
  
  t2 <- RenumberTips(t2, t1)
  t2 <- Preorder(t2)

  ME <- function(e, node) {
    expect_equal(e, 
                 .MASTSizeEdges(Postorder(t1$edge),
                                RootOnNode(t2, node = node, TRUE)$edge,
                                nTip = 5L))
  }
  
  # Values checked against brute-force rooted MAST
  ME(3L, 1)
  ME(3L, 2)
  ME(4L, 3)
  ME(3L, 4)
  ME(3L, 5)
  ME(4L, 6)
  ME(4L, 7)
  ME(3L, 8)
  ME(3L, 9)
  
  expect_equal(4L, MASTSize(t1, t2, rooted = FALSE))
})

test_that("MASTSize() matches brute force", {
  # Reported under-count: dropping t1 makes the trees identical
  tree1 <- ape::read.tree(text = "(t1,(t2,(t3,t4)),((t5,t6),(t7,t8)));")
  tree2 <- ape::read.tree(text = "(t1,(t5,t6),((t2,(t3,t4)),(t7,t8)));")
  expect_equal(MASTSize(tree1, tree2, rooted = FALSE), 7L)
  
  # Rooted: MAST lies within one child of the root of tree2
  rooted1 <- ape::read.tree(text = "((a,b),(c,(d,e)));")
  rooted2 <- ape::read.tree(text = "(e,((a,b),(c,d)));")
  expect_equal(MASTSize(rooted1, rooted2), 4L)
  expect_equal(MASTSize(rooted2, rooted1), 4L)
  
  AddRootLeaf <- function(tree) {
    ape::read.tree(text = paste0("(ROOT,", sub(";$", "", ape::write.tree(tree)),
                                 ");"))
  }
  BruteMAST <- function(tree1, tree2, rooted) {
    labels <- tree1[["tip.label"]]
    for (k in rev(seq_along(labels))[-length(labels)]) {
      for (keep in combn(labels, k, simplify = FALSE)) {
        a <- KeepTip(tree1, keep)
        b <- KeepTip(tree2, keep)
        if (rooted) {
          a <- AddRootLeaf(a)
          b <- AddRootLeaf(b)
        } else if (k < 4) {
          return(k)
        }
        if (RobinsonFoulds(a, b) == 0) return(k)
      }
    }
    1L
  }
  RandomPair <- function(nTip) {
    labels <- paste0("t", seq_len(nTip))
    Tidy <- function(tr) {
      ape::read.tree(text = ape::write.tree(RenumberTips(tr, labels)))
    }
    list(Tidy(RandomTree(labels, root = TRUE)),
         Tidy(RandomTree(labels, root = TRUE)))
  }
  
  set.seed(4110)
  for (i in 1:12) {
    pair <- RandomPair(sample(5:7, 1))
    expect_equal(MASTSize(pair[[1]], pair[[2]]),
                 BruteMAST(pair[[1]], pair[[2]], rooted = TRUE))
    expect_equal(MASTSize(pair[[1]], pair[[2]], rooted = FALSE),
                 BruteMAST(pair[[1]], pair[[2]], rooted = FALSE))
  }
})
