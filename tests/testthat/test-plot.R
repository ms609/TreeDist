library("TreeTools")

test_that("TreeDistPlot() warns", {
  expect_warning(
    expect_warning(
      expect_null(TreeDist::TreeDistPlot(PectinateTree(8))),
      "Leaves.*must be labelled with integers"),
    "fewer than 2 tips" # From plot.phylo: I don't understand why!
  )
})

test_that("TreeDistPlot() works", {
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
  skip_if(packageVersion("graphics") < "4.3")
  skip_if(packageVersion("vdiffr") < "1.0")
  
  vdiffr::expect_doppelganger("Test with space", Test1)
  vdiffr::expect_doppelganger("Test without space", Test2)
  tr$tip.label <- letters[1:11]
  vdiffr::expect_doppelganger("Test-lc-letters", Test2)
  tr$tip.label <- LETTERS[1:11]
  vdiffr::expect_doppelganger("Test-UC-LETTERS", Test2)
  
})
