library("TreeTools")

test_that("TreeDistPlot works", {
  tr <- PectinateTree(1:11)
  tr$edge.width <- rep(1:2, 10)
  Test1 <- function() {
    TreeDist::TreeDistPlot(tr, title = "Test", 
                         bold = c(2, 4, 6),
                         leaveRoom = TRUE,
                         prune = 1, graft = 10)
  }
  Test2 <- function() {
    TreeDist::TreeDistPlot(tr, title="Crop tightly", 
                           bold = c(2, 4, 6), prune = 11, graft = 10,
                           leaveRoom = FALSE)
  }
  skip_if_not_installed("vdiffr")
  skip_if(packageVersion("graphics") < "4.1")
  skip_if(packageVersion("vdiffr") < "1.0")
  
  vdiffr::expect_doppelganger("Test with space", Test1)
  vdiffr::expect_doppelganger("Test without space", Test2)
  tr$tip.label <- letters[1:11]
  vdiffr::expect_doppelganger("Test-lc-letters", Test2)
  tr$tip.label <- LETTERS[1:11]
  vdiffr::expect_doppelganger("Test-UC-LETTERS", Test2)
  
})

test_that("VisualizeMatching() works", {
  tree1 <- PectinateTree(1:11)
  tree2 <- tree1
  tree2$tip.label[c(11, 1)] <- tree1$tip.label[c(1, 11)]
  tree2r <- CollapseNode(tree2, 20:21)
  
  Minus <- function(...) {
    x <- MutualClusteringInfo(...)
    attr(x, "pairScores") <- -attr(x, "pairScores")
    x
  }
  expect_error(VisualizeMatching(Minus, PectinateTree(8), BalancedTree(8), 
                                 setPar = FALSE))
  
  skip_if_not_installed("vdiffr")
  skip_if(packageVersion("graphics") < "4.1")
  skip_if(packageVersion("vdiffr") < "1.0")
  
  TestVM <- function() {
    VisualizeMatching(MutualClusteringInfo, tree1, tree2, 
                      setPar = TRUE, precision = 3, matchZeros = FALSE,
                      Plot = plot.phylo)
  }
  vdiffr::expect_doppelganger("Test VM", TestVM)
  
  TestVMr <- function() {
    VisualizeMatching(MutualClusteringInfo, tree1, tree2r,
                      setPar = TRUE, precision = 3, matchZeros = TRUE, 
                      Plot = plot.phylo, cex = 1.5)
  }
  vdiffr::expect_doppelganger("Test VMr", TestVMr)
  
  vdiffr::expect_doppelganger("Visualize MCI matching", function() {
    par(mfrow = c(2, 2), mar = rep(0.1, 4), cex = 1.5)
    tree1 <- ape::read.tree(text="((1, 2), ((3, (4, 5)), (6, (7, (8, 9)))));")
    tree2 <- ape::read.tree(text="((1, 2), ((3, 4, (5, 9)), (6, (7, 8))));")
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
  
  vdiffr::expect_doppelganger("RF: Collapse a node", function() {
    par(mfrow = c(2, 2), mar = rep(0.1, 4), cex = 1.5)
    tree1 <- ape::read.tree(text="((1, 2), ((6, (7, 8)), (3, 4, (5, 9))));")
    tree2 <- ape::read.tree(text="((1, 2), ((3, (4, (5, 9))), (6, (7, 8))));")
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
  
  
  vdiffr::expect_doppelganger("RF: Collapse and change", function() {
    par(mfrow = c(2, 2), mar = rep(0.1, 4), cex = 1.5)
    tree1 <- ape::read.tree(text="((1, 2), ((6, (7, 8)), (3, 4, (5, 9))));")
    tree2 <- ape::read.tree(text="((1, 2), ((3, (4, (5, 9))), ((6, 7), 8)));")
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
  
  vdiffr::expect_doppelganger("RF VM Single splits; plainEdges", function() {
    par(mfrow = c(2, 2), mar = rep(0.1, 4), cex = 1.5)
    tree1 <- ape::read.tree(text = "((1, 2), (3, 4, 5, 6, 7, 8));")
    tree2 <- ape::read.tree(text = "((1, 2, 3), (4, 5, 6, 7, 8));")
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
  
  vdiffr::expect_doppelganger("JRF VM matchZeros FALSE", function() {
    JRF2 <- function(tree1, tree2, ...) 
      JaccardRobinsonFoulds(tree1, tree2, k = 2, allowConflict = FALSE, ...)
    
    tree1 <- EnforceOutgroup(as.phylo(704564, 10), paste0("t", c(1,4,5,8,9)))
    tree2 <- EnforceOutgroup(as.phylo(20165 , 10), paste0("t", c(1,4)))
    VisualizeMatching(JRF2, tree1, tree2, matchZeros = FALSE)
  })
})

test_that("MST example plots as expected", {
  skip_if_not_installed("graphics", "4.1")
  skip_if_not_installed("vdiffr", "1.0")
  skip_if_not_installed("TreeTools", "1.6.0.9008")
  vdiffr::expect_doppelganger("MST example plot", function() {
    distances <- structure(
      c(3, 2.3, 2.3, 2.3, 3, 1.7, 2.3, 2.3, 2.3, 4.5, 4.5, 1.7, 1.7, 2.3, 3,
        2.3, 1.7, 2.3, 2.3, 3.8, 4.4, 1.6, 3.1, 2.3, 2.6, 3, 1.5, 2.6, 4.7,
        4.6, 2.6, 2.3, 1.5, 3, 2.6, 1.5, 4.3, 4.7, 1.7, 2.6, 1.5, 3, 1.6, 4.6,
        4.7, 2.3, 2.3, 1.7, 1.7, 4.4, 3.8, 3.1, 3.1, 1.5, 4.4, 4.4, 3.1, 2.6,
        4.2, 4.8, 3, 4.8, 4.2, 4.7, 4.3, 1.6),
      class = "dist", Size = 12L, Diag = FALSE, Upper = FALSE
    ) # dput(round(ClusteringInfoDist(as.phylo(5:16, 8)), 1))
    
    mapping <- structure(
      c(0.75, 0.36, 0.89, 0.85, 0.89, 0.36, 0.64, 0.6, 0.6, 0.85, -3.4, -3.4,
        0, -0.25, 1.33, 0.39, -1.33, 0.25, 0, -1.52, 1.52, -0.39, -0.58, 0.58),
      dim = c(12L, 2L)
    ) # dput(round(cmdscale(distances, k = 2), 2))
    
    mstEnds <- MSTEdges(distances)
    
    # Set up blank plot
    plot(mapping, asp = 1, frame.plot = FALSE, ann = FALSE, axes = FALSE,
         type = "n")
    
    # Add MST
    MSTSegments(mapping, mstEnds, 
                col = StrainCol(distances, mapping, mstEnds))
    
    # Add points at end so they overprint the MST
    points(mapping)
    SpectrumLegend(legend = c("Contracted", "Median", "Extended"),
                   palette = hcl.colors(256L, "RdYlBu", rev = TRUE))
  })
})

test_that("StrainCol() handles zeroes", {
  distances <- dist(c(1, 1, 10, 100))
  mapping <- cmdscale(distances)
  expect_equal(attr(StrainCol(distances, mapping), "logStrain")[1], Inf)
})
