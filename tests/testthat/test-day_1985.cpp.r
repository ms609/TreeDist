context("Day 1985")

test_that("Day 1985 examples", {
  library("TreeTools", quietly = TRUE, warn.conflicts = FALSE)
  
  PrepareTree <- function (text) {
    tmp <- ape::read.tree(text = text)
    #Preorder(RenumberTips(RootTree(tmp, '1'), as.character(seq_along(tmp$tip.label))))
    Preorder(RenumberTips(tmp, as.character(seq_along(tmp$tip.label))))
  }
  t1 <- PrepareTree("((10, 7), (6, (8, 11)), (12, (4, (2, 1))), 14, (5, 9, 13), 3);")
  plot(t1); edgelabels(); nodelabels()
  t2 <- PrepareTree("(((2, 4, 5, 7, 9, 10, 12, 13), (1, 14)), (6, (8, 11)), 3);")
  
  COMCLUST(list(t1, t2))
  expect_null(NULL)
  t1n <- structure(list(edge = structure(
    c(22, 15, 15, 22, 19, 19, 16, 16, 22, 21, 21, 20, 20, 17, 17, 22, 22, 18, 18, 
      18, 22, 15, 10L, 7L, 19, 6L, 16, 8L, 11L, 21, 12L, 20, 
      4L, 17, 2L, 1L, 14L, 18, 5L, 9L, 13L, 3L), .Dim = c(21L, 2L)),
    Nnode = 8L, tip.label = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14")), class = "phylo")
  
})