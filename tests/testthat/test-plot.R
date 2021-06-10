library('TreeTools')
context("plot.R")

if (packageVersion('graphics') < '4.1')
test_that('TreeDistPlot works', {
  tr <- PectinateTree(1:11)
  tr$edge.width <- rep(1:2, 10)
  Test1 <- function () {
    TreeDist::TreeDistPlot(tr, title='Test', 
                         bold = c(2, 4, 6),
                         leaveRoom = TRUE,
                         prune = 1, graft = 10)
  }
  Test2 <- function () {
    TreeDist::TreeDistPlot(tr, title='Crop tightly', 
                           bold = c(2, 4, 6), prune = 11, graft = 10,
                           leaveRoom = FALSE)
  }
  skip_if_not_installed('vdiffr')
  skip_if(packageVersion("graphics") < "4.1")
  skip_if(packageVersion("vdiffr") < "1.0")
  
  vdiffr::expect_doppelganger("Test with space", Test1)
  vdiffr::expect_doppelganger("Test without space", Test2)
  tr$tip.label <- letters[1:11]
  vdiffr::expect_doppelganger("Test-lc-letters", Test2)
  tr$tip.label <- LETTERS[1:11]
  vdiffr::expect_doppelganger("Test-UC-LETTERS", Test2)
  
}) else expect_false(packageVersion('graphics') < '4.1')

if (packageVersion('graphics') < '4.1')
test_that('VisualizeMatching() works', {
  tree1 <- PectinateTree(1:11)
  tree2 <- tree1
  tree2$tip.label[c(11, 1)] <- tree1$tip.label[c(1, 11)]
  tree2r <- CollapseNode(tree2, 20:21)
  
  Minus <- function (...) {
    x <- MutualClusteringInfo(...)
    attr(x, 'pairScores') <- -attr(x, 'pairScores')
    x
  }
  expect_error(VisualizeMatching(Minus, PectinateTree(8), BalancedTree(8), 
                                 setPar = FALSE))
  
  skip_if_not_installed('vdiffr')
  skip_if(packageVersion("graphics") < "4.1")
  skip_if(packageVersion("vdiffr") < "1.0")
  
  TestVM <- function () {
    VisualizeMatching(MutualClusteringInfo, tree1, tree2, 
                      setPar = TRUE, precision = 3, matchZeros = FALSE,
                      Plot = plot.phylo)
  }
  vdiffr::expect_doppelganger('Test VM', TestVM)
  
  TestVMr <- function () {
    VisualizeMatching(MutualClusteringInfo, tree1, tree2r,
                      setPar = TRUE, precision = 3, matchZeros = TRUE, 
                      Plot = plot.phylo, cex = 1.5)
  }
  vdiffr::expect_doppelganger('Test VMr', TestVMr) # Unclear why this test fails on Travis. 
  
  vdiffr::expect_doppelganger('Visualize MCI matching', function () {
    par(mfrow = c(2, 2), mar = rep(0.1, 4), cex = 1.5)
    tree1 <- ape::read.tree(text='((1, 2), ((3, (4, 5)), (6, (7, (8, 9)))));')
    tree2 <- ape::read.tree(text='((1, 2), ((3, 4, (5, 9)), (6, (7, 8))));')
    VisualizeMatching(MutualClusteringInfo, tree1, tree2,
                      setPar = FALSE, precision = 3L,
                      Plot = TreeDistPlot,
                      matchZeros = FALSE,
                      leaveRoom = FALSE)
    VisualizeMatching(MutualClusteringInfo, tree2, tree1,
                      setPar = FALSE, precision = 3,
                      Plot = TreeDistPlot,
                      matchZeros = FALSE,
                      leaveRoom = FALSE)
  })
  
  vdiffr::expect_doppelganger('RF: Collapse a node', function () {
    par(mfrow = c(2, 2), mar = rep(0.1, 4), cex = 1.5)
    tree1 <- ape::read.tree(text='((1, 2), ((6, (7, 8)), (3, 4, (5, 9))));')
    tree2 <- ape::read.tree(text='((1, 2), ((3, (4, (5, 9))), (6, (7, 8))));')
    VisualizeMatching(RobinsonFouldsMatching, tree1, tree2,
                      setPar = FALSE, precision = 3,
                      Plot = TreeDistPlot,
                      matchZeros = TRUE,
                      leaveRoom = FALSE)
    VisualizeMatching(RobinsonFouldsMatching, tree2, tree1,
                      setPar = FALSE, precision = 3,
                      Plot = TreeDistPlot,
                      matchZeros = FALSE,
                      leaveRoom = FALSE)
  })
  
  
  vdiffr::expect_doppelganger('RF: Collapse and change', function () {
    par(mfrow = c(2, 2), mar = rep(0.1, 4), cex = 1.5)
    tree1 <- ape::read.tree(text='((1, 2), ((6, (7, 8)), (3, 4, (5, 9))));')
    tree2 <- ape::read.tree(text='((1, 2), ((3, (4, (5, 9))), ((6, 7), 8)));')
    VisualizeMatching(RobinsonFouldsMatching, tree1, tree2,
                      setPar = FALSE, precision = 3L,
                      Plot = TreeDistPlot,
                      matchZeros = TRUE,
                      leaveRoom = FALSE)
    VisualizeMatching(RobinsonFouldsMatching, tree2, tree1,
                      setPar = FALSE, precision = 3L,
                      Plot = TreeDistPlot,
                      matchZeros = FALSE,
                      leaveRoom = FALSE)
  })
  
  vdiffr::expect_doppelganger('RF VM Single splits; plainEdges', function () {
    par(mfrow = c(2, 2), mar = rep(0.1, 4), cex = 1.5)
    tree1 <- ape::read.tree(text = '((1, 2), (3, 4, 5, 6, 7, 8));')
    tree2 <- ape::read.tree(text = '((1, 2, 3), (4, 5, 6, 7, 8));')
    VisualizeMatching(RobinsonFouldsMatching, tree1, tree2,
                      setPar = FALSE,
                      Plot = TreeDistPlot,
                      matchZeros = TRUE,
                      plainEdges = TRUE,
                      edge.width = NULL,
                      leaveRoom = FALSE)
    VisualizeMatching(RobinsonFouldsMatching, tree2, tree1,
                      setPar = FALSE,
                      Plot = TreeDistPlot,
                      matchZeros = FALSE,
                      plainEdges = FALSE,
                      leaveRoom = FALSE)
  })
  
  vdiffr::expect_doppelganger('JRF VM matchZeros FALSE', function () {
    JRF2 <- function (tree1, tree2, ...) 
      JaccardRobinsonFoulds(tree1, tree2, k = 2, allowConflict = FALSE, ...)
    
    tree1 <- EnforceOutgroup(as.phylo(704564, 10), paste0('t', c(1,4,5,8,9)))
    tree2 <- EnforceOutgroup(as.phylo(20165 , 10), paste0('t', c(1,4)))
    VisualizeMatching(JRF2, tree1, tree2, matchZeros = FALSE)
  })
})