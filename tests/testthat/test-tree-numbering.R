library(ape)

#' @ author Martin R. Smith <martins@gmail.com>
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

context("Test tree rearrangement")
small_tree <- rtree(8)
large_tree <- rtree(80)  
test_that("NNI trees conform to phylo expectations", {
  for (i in 1:100)  CheckTreeSanity(small_tree <- NNI(small_tree))
  for (i in 1:1000) CheckTreeSanity(large_tree <- NNI(large_tree))
})
test_that("SPR trees conform to phylo expectations", {
  for (i in 1:100)  CheckTreeSanity(small_tree <- SPR(small_tree))
  for (i in 1:1000) CheckTreeSanity(large_tree <- SPR(large_tree))
})
test_that("TBR trees conform to phylo expectations", {
  for (i in 1:100)  CheckTreeSanity(small_tree <- TBR(small_tree))
  for (i in 1:1000) CheckTreeSanity(large_tree <- TBR(large_tree))
})
