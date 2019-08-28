library('vdiffr')
# use manage_cases() to manage
context("plot.R")

test_that('TreeDistPlot works', {
  tr <- ape::read.tree(text='(1, (2, (3, (4, (5, (6, (7, (8, (9, (10, 11))))))))));')
  Test1 <- function () {
    TreeDist::TreeDistPlot(tr, title='Test', 
                         bold=c(2, 4, 6),
                         leaveRoom = TRUE,
                         prune=1, graft=10)
  }
  Test2 <- function () {
    TreeDist::TreeDistPlot(tr, title='Crop tightly', 
                                              bold=c(2, 4, 6),
                                              leaveRoom = FALSE,
                                              prune=11, graft=10)
  }
  expect_doppelganger("Test with space", Test1)
  expect_doppelganger("Test without space", Test2)
})

test_that('VisualizeMatching works', {
  tree1 <-  ape::read.tree(text='(1, (2, (3, (4, (5, (6, (7, (8, (9, (10, 11))))))))));')
  tree2 <-  ape::read.tree(text='(11, (2, (3, (4, (5, (6, (7, (8, (9, (10, 1))))))))));')
  TestVM <- function () {
    VisualizeMatching(MutualClusteringInfo, tree1, tree2, 
                                          setPar = TRUE, precision=3, 
                                          Plot = plot.phylo)
  }
  expect_doppelganger('Test VM', TestVM)
})