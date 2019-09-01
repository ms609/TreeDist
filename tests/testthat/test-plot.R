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
  tree2r <-  ape::read.tree(text='(11, (2, (3, (4, (5, (6, (7, (8, 9, 10, 1))))))));')
  TestVM <- function () {
    VisualizeMatching(MutualClusteringInfo, tree1, tree2, 
                                          setPar = TRUE, precision=3, 
                                          Plot = plot.phylo)
  }
  TestVMr <- function () {
  }
  expect_doppelganger('Test VM', TestVM)
  expect_doppelganger('Test VMr', TestVMr)
  expect_doppleganger('Collapse a node', function () {
    par(mfrow=c(2, 2))
    tree1 <- ape::read.tree(text='((1, 2), ((6, (7, 8)), (3, 4, (5, 9))));')
    tree2 <- ape::read.tree(text='((1, 2), ((3, (4, (5, 9))), (6, (7, 8))));')
    VisualizeMatching(RobinsonFoulds, tree1, tree2,
                      setPar = FALSE, precision=3,
                      Plot = TreeDistPlot,
                      leaveRoom=FALSE)
    VisualizeMatching(RobinsonFoulds, tree2, tree1,
                      setPar = FALSE, precision=3,
                      Plot = TreeDistPlot,
                      leaveRoom=FALSE)
  })
  
  
  expect_doppleganger('Collapse a node', function () {
    par(mfrow=c(2, 2))
    tree1 <- ape::read.tree(text='((1, 2), ((6, (7, 8)), (3, 4, (5, 9))));')
    tree2 <- ape::read.tree(text='((1, 2), ((3, (4, (5, 9))), ((6, 7), 8)));')
    VisualizeMatching(RobinsonFoulds, tree1, tree2,
                      setPar = FALSE, precision=3,
                      Plot = TreeDistPlot,
                      leaveRoom=FALSE)
    VisualizeMatching(RobinsonFoulds, tree2, tree1,
                      setPar = FALSE, precision=3,
                      Plot = TreeDistPlot,
                      leaveRoom=FALSE)
  })
  
  
})