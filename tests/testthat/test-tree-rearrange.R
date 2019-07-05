library(ape)

context("Tree rearrangements")
tree5a <- read.tree(text='(a, (b, (c, (d, e))));')
tree5b <- read.tree(text='((a, b), (c, (d, e)));')
tree6  <- Preorder(read.tree(text="((a, (b, (c, d))), (e, f));"))
tree6b <- Preorder(read.tree(text="((a, (b, c)), (d, (e, f)));"))
tree8  <- read.tree(text="(((a, (b, (c, d))), (e, f)), (g, h));")
tree11 <- read.tree(text="((((a, b), (c, d)), e), ((f, (g, (h, i))), (j, k)));")
attr(tree5a, 'order') <- attr(tree5b, 'order') <- attr(tree8, 'order') <- attr(tree11, 'order') <- 'preorder'

test_that("Malformed trees throw errors", {
  treeDoubleNode    <- read.tree(text = "((((((1,2)),3),4),5),6);")
  treePolytomy      <- read.tree(text = "((((1,2,3),4),5),6);")
  treeDoublyPoly    <- read.tree(text = "(((((1,2,3)),4),5),6);")
  
  expect_error(Preorder(treeDoubleNode))
  expect_error(Preorder(treePolytomy))
  expect_error(Preorder(treeDoublyPoly))
  
  expect_error(NNI(treeDoubleNode))
  expect_error(NNI(treePolytomy))
  expect_error(NNI(treeDoublyPoly))
  
  expect_error(SPR(treeDoubleNode))
  expect_error(SPR(treePolytomy))
  expect_error(SPR(treeDoublyPoly))
  
  expect_error(TBR(treeDoubleNode))
  expect_error(TBR(treePolytomy))
  expect_error(TBR(treeDoublyPoly))
  
})

test_that("NNI works", {
  trComb    <- read.tree(text = "(((((1,2),3),4),5),6);")
  suppressWarnings(RNGversion("3.5.0")) # Until we can require R3.6.0
  set.seed(0)
  nniComb <- NNI(trComb)
  expect_equal(nniComb$tip.label, trComb$tip.label)
  expect_equal(nniComb$Nnode, trComb$Nnode)
  expect_equal(nniComb, read.tree(text = "(((((3,2),1),4),5),6);"))  
})

test_that("TBR can swap over root", {
  expect_equal(TBR(tree5a, 1, c(7, 1)), read.tree(text='(a, (d, (e, (c, b))));'))
  expect_equal(TBR(tree5a, 2, c(5, 1)), read.tree(text='(a, (c, (b, (d, e))));'))
  expect_equal(TBR(tree5b, 1, c(7, 1)), read.tree(text='((a, b), (d, (c, e)));'))
  expect_equal(TBR(tree5b, 4, c(7, 1)), read.tree(text='((a, b), (d, (c, e)));'))
})

test_that("TBR works", {
  tree <- tree8
  ### expect_equal(TBR(tree, 3, 1 ), read.tree(text="((a, ((b, (c, d)), (e, f))), (g, h));"))
  ### expect_warning(expect_identical(TBR(tree, 3, 2), tree))
  ### expect_warning(expect_identical(TBR(tree, 3, 3), tree))
  ### expect_warning(expect_identical(TBR(tree, 3, 4), tree))
  ### expect_warning(expect_identical(TBR(tree, 3, 44), tree))
  ### expect_equal(TBR(tree, 3, 5 ), read.tree(text="((((a, b), (c, d)), (e, f)), (g, h));"))
  ### expect_equal(TBR(tree, 3, 6 ), read.tree(text="(((b, (a, (c, d))), (e, f)), (g, h));"))
  ### expect_equal(TBR(tree, 3, 7 ), read.tree(text="(((b, ((a, c), d)), (e, f)), (g, h));"))
  ### expect_equal(TBR(tree, 3, 8 ), read.tree(text="(((b, (c, (a, d))), (e, f)), (g, h));"))
  ### expect_equal(TBR(tree, 3, 9 ), read.tree(text="(((b, (c, d)), (a, (e, f))), (g, h));"))
  ### expect_equal(TBR(tree, 3, 10), read.tree(text="(((b, (c, d)), ((a, e), f)), (g, h));"))
  ### expect_equal(TBR(tree, 3, 11), read.tree(text="(((b, (c, d)), (e, (a, f))), (g, h));"))
  ### expect_equal(TBR(tree, 3, 12), read.tree(text="(((b, (c, d)), (e, f)), (a, (g, h)));"))
  ### expect_equal(TBR(tree, 3, 13), read.tree(text="(((b, (c, d)), (e, f)), ((g, a), h));"))
  ### expect_equal(TBR(tree, 3, 14), read.tree(text="(((b, (c, d)), (e, f)), (g, (a, h)));"))
  
  tree <- tree8
  expect_equal(TBR(tree, 6, c(1 , 6)), read.tree(text="((((a, b), (e, f)), (c, d)), (g, h));"))
  expect_equal(TBR(tree, 6, c(1 , 7)), read.tree(text="((((a, b), (e, f)), (c, d)), (g, h));"))
  expect_equal(TBR(tree, 6, c(1 , 8)), read.tree(text="((((a, b), (e, f)), (c, d)), (g, h));"))
  expect_equal(TBR(tree, 6, c(2 , 6)), TBR(tree, 6, c(2 , 7)))
  expect_equal(TBR(tree, 6, c(2 , 6)), TBR(tree, 6, c(2 , 8)))
  expect_equal(TBR(tree, 6, c(2 , 6)), read.tree(text="((((a, b), (c, d)), (e, f)), (g, h));"))
  expect_equal(TBR(tree, 6, c(3 , 6)), read.tree(text="(((((c, d), a), b), (e, f)), (g, h));"))
  expect_warning(expect_identical(TBR(tree, 6, c(4 , 6)), tree))
  expect_warning(expect_identical(TBR(tree, 8, c(6 , 8)), tree))
  expect_warning(expect_identical(TBR(tree, 6, c(5 , 6)), tree))
  expect_warning(expect_identical(TBR(tree, 6, c(6 , 6)), tree))
  expect_warning(expect_identical(TBR(tree, 6, c(6 , 7)), tree))
  expect_warning(expect_identical(TBR(tree, 6, c(6 , 8)), tree))
  expect_equal(TBR(tree, 6, c(9 , 6)), read.tree(text="(((a, b), ((c, d), (e, f))), (g, h));"))
  expect_equal(TBR(tree, 6, c(10, 6)), read.tree(text="(((a, b), (((c, d), e), f)), (g, h));"))
  expect_equal(TBR(tree, 6, c(11, 6)), read.tree(text="(((a, b), (((c, d), f), e)), (g, h));"))
  expect_equal(TBR(tree, 6, c(12, 6)), read.tree(text="(((a, b), (e, f)), ((c, d), (g, h)));"))
  expect_equal(TBR(tree, 6, c(13, 6)), read.tree(text="(((a, b), (e, f)), (((c, d), g), h));"))
  expect_equal(TBR(tree, 6, c(14, 6)), read.tree(text="(((a, b), (e, f)), (((c, d), h), g));"))
  expect_warning(expect_identical(TBR(tree, 6, c(6, 15)), tree))
  
  expect_equal(TBR(tree, 4, c(1, 5)),  read.tree(text="(((a, (e, f)), (b, (c, d))), (g, h));"))
  expect_equal(TBR(tree, 4, c(1, 6)),  read.tree(text="(((a, (e, f)), (b, (c, d))), (g, h));"))
  expect_equal(TBR(tree, 4, c(1, 7)),  read.tree(text="(((a, (e, f)), (c, (b, d))), (g, h));"))
  expect_equal(TBR(tree, 4, c(1, 8)),  read.tree(text="(((a, (e, f)), (d, (b, c))), (g, h));"))
  
  tree <- tree11 
  tree$edge.length = rep(1, 20) 
  expect_equal(TBR(tree11, 11, c(8, 17)), read.tree(text='((j, k), (e, ((a, b), (c, (d, (i, (h, (g, f))))))));'))
  expect_equal(TBR(tree11, 11, c(2, 11)), read.tree(text='((j, k), (e, (((a, b), (c, d)), (f, (g, (i, h))))));'))
  expect_warning(TBR(tree11, 10, c(2, 11)))
  expect_equal(TBR(tree11, 10, c(3, 11)), read.tree(text='(e, ((c, d), ((a, b), ((j, k), (f, (g, (h, i)))))));'))
    
})

test_that("RootedTBR fails", {
  #  tree8 <- read.tree(text="(((a, (b, (c, d))), (e, f)), (g, h));")
  #  tree11 <- read.tree(text="((((a, b), (c, d)), e), ((f, (g, (h, i))), (j, k)));")

  expect_equal(TBR(tree8, 4, c(3, 7)), RootedTBR(tree8, 4, c(3, 7)))
  expect_equal(TBR(tree8, 4, c(1, 5)), RootedTBR(tree8, 4, c(1, 5)))
  expect_warning(RootedTBR(tree5a, edgeToBreak = 1))
  expect_warning(RootedTBR(tree5a, edgeToBreak = 2))
  expect_equal(RootedTBR(tree5a, edgeToBreak = 3, mergeEdges=6), read.tree(text='(a, (c, (b, (d, e))));'))
  expect_silent(replicate(100, RootedTBR(tree5a)))
  expect_warning(RootedTBR(tree8, 4, c(13, 6)))
  expect_warning(RootedTBR(read.tree(text='((a, b), (c, d));')))
})

test_that("RootedSPR fails", {
  expect_warning(RootedSPR(read.tree(text='((a, b), (c, d));')))
  expect_warning(RootedSPR(tree8, edgeToBreak=1))
  expect_warning(RootedSPR(tree8, edgeToBreak=13))
  expect_warning(RootedSPR(tree8, edgeToBreak=14))
  warnTree1 <- read.tree(text='((a, (b, (c, d))), (e, (f, (g, h))));')
  warnTree2 <- read.tree(text='((a, (b, (c, d))), (((e, f), g), h));')
  attr(warnTree1, 'order') <- attr(warnTree2, 'order') <- 'preorder'
  expect_warning(RootedSPR(warnTree1, 3))
  expect_warning(RootedSPR(warnTree1, 10))
  expect_warning(RootedSPR(warnTree2, 9))
  expect_warning(RootedSPR(warnTree2, 8))
})

test_that("SPR is special case of TBR", {
  #library(devtools); library(testthat); library(ape); load_all()
  Plot <- function (x) {plot(x); nodelabels(cex=0.8); edgelabels()}
  
  expect_equal(SPR(tree11, 3, 9), TBR(tree11, 3, c(3, 9)))
  expect_equal(SPR(tree11, 12, 9), TBR(tree11, 12, c(12, 9)))
  expect_equal(root(SPR(tree11, 1, 14), letters[1:5], resolve.root=TRUE), TBR(tree11, 1, c(1, 14)))
  expect_error(SPR(tree11, 1, 6))
})

test_that("TBR move lister works", {
  edge <- tree6$edge
  parent <- edge[, 1]
  child <- edge[, 2]
  moves <- TBRMoves(parent, child)
  expect_equal(rep(2:10, c(7, 5, 6, 4, 6, 6, 4, 6, 6)), moves[, 1])
  expect_equal(c(4:10, 6:10, 2, 6:10, 2, 8:10, 
                 rep(c(2:4, 8:10), 2), 4:7, 2:7, 2:7), moves[, 2])
  rootedMoves <- TBRMoves(parent, child, retainRoot=TRUE)
  expect_equal(matrix(c(2,4, 2,5, 2,6, 2,7,
                        3,6, 3,7,
                        4,2, 4,6, 4,7,
                        5,2,
                        6,2, 6,3, 6,4, 
                        7,2, 7,3, 7,4), ncol=2, byrow=TRUE), rootedMoves)
  
  edge <- tree6b$edge
  parent <- edge[, 1]
  child <- edge[, 2]
  rootedMoves <- TBRMoves(parent, child, retainRoot=TRUE)
  expect_equal(matrix(c(2,2,4,5,7,7 ,9,10,
                        4,5,2,2,9,10,7,7 ), ncol=2), rootedMoves)
  expect_equal(length(AllTBR(parent, child, retainRoot=TRUE)), 4)

})

test_that("CollapseNodes works", {
  
  expect_error(CollapseNode(1:5, tree8))
  expect_error(CollapseNode(tree8, 1))
  
  suppressWarnings(RNGversion("3.5.0")) # Until we can require R3.6.0
  set.seed(1)
  
  tree <- rtree(7)
  expect_equal(tree, CollapseNode(tree, integer(0)))
  
  no1213 <- CollapseNode(tree, c(12, 13))
  expect_equal(no1213$edge, matrix(c(8, 9, 9, 8, 10, 11, 11, 10, 10, 10, 
                                      9, 1, 2, 10, 11, 3:7), ncol=2))
  el <- tree$edge.length
  expect_equal(no1213$edge.length, c(el[1:7], el[8] + c(c(el[9] + el[10:11]), el[12])))
  
  no11 <- CollapseEdge(tree, 5L)
  expect_equal(no11$edge, matrix(c(8, 9, 9, 8, 10, 10, 10, 11, 12, 12, 11,
                                   9, 1, 2, 10, 3, 4, 11, 12, 5:7), ncol=2))

})

#' @author Martin R. Smith <martins@gmail.com>
CheckTreeSanity <- function (tree) {
  nTip <- length(tree$tip.label)
  nNode <- tree$Nnode
  edge <- tree$edge
  parent <- edge[, 1]
  child <- edge[, 2]
  aok <- TRUE
  expect_true(all(parent > nTip),
              info=paste0("Parent nodes on edge(s) ", paste(which(parent <= nTip), collapse=', '), 
                          " are tips (nTip = ", nTip, ')')
  )
  expect_equal(min(parent), nTip + 1,
               info=paste0("Root is numbered ", min(parent), "; expecting ", nTip + 1)
  )
  expect_false(min(parent) %in% child, 
               info=paste0("Root node (", min(parent), ") is child of edge ", paste0(which(min(parent) == child), collapse=', '))
  )
  expect_true(all(seq_len(nTip) %in% child)) # No missing tips
  expect_equal(max(parent), nTip + nNode)
  tips <- child <= nTip
  expect_equal(sum(tips), nTip)
  expect_true(all(child[!tips] > parent[!tips]), info="Parent nodes must be > child nodes")
}

suppressWarnings(RNGversion("3.5.0")) # Until we can require R3.6.0
set.seed(0)
small_tree <- rtree(8)
large_tree <- rtree(80)  
test_that("NNI trees conform to phylo expectations", {
  for (i in 1:60)  CheckTreeSanity(small_tree <- NNI(small_tree))
  for (i in 1:250) CheckTreeSanity(large_tree <- NNI(large_tree))
})
test_that("SPR trees conform to phylo expectations", {
  for (i in 1:60)  CheckTreeSanity(small_tree <- SPR(small_tree))
  for (i in 1:250) CheckTreeSanity(large_tree <- SPR(large_tree))
})
test_that("TBR trees conform to phylo expectations", {
  for (i in 1:60)  CheckTreeSanity(small_tree <- TBR(small_tree))
  for (i in 1:250) CheckTreeSanity(large_tree <- TBR(large_tree))
})
