library(ape)

context("tree rearrangement")
test_that("NNI works", {
  trComb    <- read.tree(text = "(((((1,2),3),4),5),6);")
  set.seed(0)
  nniComb <- NNI(trComb)
  expect_equal(nniComb$tip.label, trComb$tip.label)
  expect_equal(nniComb$Nnode, trComb$Nnode)
  
  
  trCombLen <- read.tree(text = "(((((1:1,2:1):1,3:2):1,4:3):1,5:1):1,6:1);")

  expect_equal(NNI(trComb), read.tree(text = "(((((3,2),1),4),5),6);"))
  
})
context("tree rearrangement")
test_that("TBR works", {
  tree <- read.tree(text="(((a, (b, (c, d))), (e, f)), (g, h));")
  expect_equal(TBR(tree, 3, 1 ),read.tree(text="((a, ((b, (c, d)), (e, f))), (g, h));"))
  expect_warning(expect_identical(TBR(tree, 3, 2), tree))
  expect_warning(expect_identical(TBR(tree, 3, 3), tree))
  expect_warning(expect_identical(TBR(tree, 3, 4), tree))
  expect_equal(TBR(tree, 3, 5 ), read.tree(text="((((a, b), (c, d)), (e, f)), (g, h));"))
  expect_equal(TBR(tree, 3, 6 ), read.tree(text="(((b, (a, (c, d))), (e, f)), (g, h));"))
  expect_equal(TBR(tree, 3, 7 ), read.tree(text="(((b, ((a, c), d)), (e, f)), (g, h));"))
  expect_equal(TBR(tree, 3, 8 ), read.tree(text="(((b, (c, (a, d))), (e, f)), (g, h));"))
  expect_equal(TBR(tree, 3, 9 ), read.tree(text="(((b, (c, d)), (a, (e, f))), (g, h));"))
  expect_equal(TBR(tree, 3, 10), read.tree(text="(((b, (c, d)), ((a, e), f)), (g, h));"))
  expect_equal(TBR(tree, 3, 11), read.tree(text="(((b, (c, d)), (e, (a, f))), (g, h));"))
  expect_equal(TBR(tree, 3, 12), read.tree(text="(((b, (c, d)), (e, f)), (a, (g, h)));"))
  expect_equal(TBR(tree, 3, 13), read.tree(text="(((b, (c, d)), (e, f)), ((g, a), h));"))
  expect_equal(TBR(tree, 3, 14), read.tree(text="(((b, (c, d)), (e, f)), (g, (a, h)));"))
  
  
  expect_equal(TBR(tree, 6, c(1 , 6)), read.tree(text="((((a, b), (e, f)), (c, d)), (g, h));"))
  expect_equal(TBR(tree, 6, c(1 , 7)), read.tree(text="((((a, b), (e, f)), (c, d)), (g, h));"))
  expect_equal(TBR(tree, 6, c(1 , 8)), read.tree(text="((((a, b), (e, f)), (c, d)), (g, h));"))
  expect_equal(TBR(tree, 6, c(2 , 6)), read.tree(text="(((a, b), (e, f)), (g, h));"))
  expect_equal(TBR(tree, 6, c(3 , 6)), read.tree(text="(((a, b), (e, f)), (g, h));"))
  expect_equal(TBR(tree, 6, c(4 , 6)), read.tree(text="(((a, b), (e, f)), (g, h));"))
  expect_equal(TBR(tree, 6, c(5 , 6)), read.tree(text="(((a, b), (e, f)), (g, h));"))
  expect_equal(TBR(tree, 6, c(6 , 6)), read.tree(text="(((a, b), (e, f)), (g, h));"))
  expect_equal(TBR(tree, 6, c(7 , 6)), read.tree(text="(((a, b), (e, f)), (g, h));"))
  expect_equal(TBR(tree, 6, c(8 , 6)), read.tree(text="(((a, b), (e, f)), (g, h));"))
  expect_equal(TBR(tree, 6, c(9 , 6)), read.tree(text="(((a, b), (e, f)), (g, h));"))
  expect_equal(TBR(tree, 6, c(10, 6)), read.tree(text="(((a, b), (e, f)), (g, h));"))
  expect_equal(TBR(tree, 6, c(11, 6)), read.tree(text="(((a, b), (e, f)), (g, h));"))
  expect_equal(TBR(tree, 6, c(12, 6)), read.tree(text="(((a, b), (e, f)), (g, h));"))
  expect_equal(TBR(tree, 6, c(13, 6)), read.tree(text="(((a, b), (e, f)), (g, h));"))
  expect_equal(TBR(tree, 6, c(14, 6)), read.tree(text="(((a, b), (e, f)), (g, h));"))
  
  expect_equal(TBR(tree, 4, c(1, 5)),  read.tree(text="(((a, (e, f)), (b, (c, d))), (g, h));"))
  expect_equal(TBR(tree, 4, c(1, 6)),  read.tree(text="(((a, (e, f)), (b, (c, d))), (g, h));"))
  expect_equal(TBR(tree, 4, c(1, 7)),  read.tree(text="(((a, (e, f)), (c, (b, d))), (g, h));"))
  expect_equal(TBR(tree, 4, c(1, 8)),  read.tree(text="(((a, (e, f)), (d, (b, c))), (g, h));"))

                                                
  
  TBRit(tree, 4, c(1 , 7))
  
  
})
