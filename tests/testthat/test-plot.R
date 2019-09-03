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
  skip_on_travis() # Skips all following tests in this block
  expect_doppelganger('Test VMr', TestVMr) # Unclear why this test fails on Travis. 
  expect_doppelganger('Collapse a node', function () {
    par(mfrow=c(2, 2), mar = rep(0.1, 4), cex=1.5)
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
  
  
  expect_doppelganger('Collapse and change', function () {
    par(mfrow=c(2, 2), mar = rep(0.1, 4), cex=1.5)
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
  
  expect_doppelganger('VM Single splits', function () {
    par(mfrow=c(2, 2), mar = rep(0.1, 4), cex=1.5)
    tree1 <- ape::read.tree(text='((1, 2), (3, 4, 5, 6, 7, 8));')
    tree2 <- ape::read.tree(text='((1, 2, 3), (4, 5, 6, 7, 8));')
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