library("TreeTools")

test_that("as.phylo.HPart", {
  bal7 <- BalancedTree(7)
  hb7 <- as.HPart(bal7)
  expect_equal(Preorder(as.phylo.HPart(hb7)), bal7)
  
  bal17 <- BalancedTree(17)
  hb17 <- as.HPart(bal17)
  expect_equal(Preorder(as.phylo.HPart(hb17)), bal17)
})

test_that("HParts are relabelled correctly", {
  bal7 <- BalancedTree(7)
  hb7 <- as.HPart(bal7)
  
  map <- c(7:4, 1:3)
  mappedLabels <- paste0("t", map)
  
  hbMap <- RenumberTips(hb7, mappedLabels)
  # Here we want only to map the internal node IDs
  attr(hbMap, "tip.label") <- TipLabels(hb7)
  
  bal7tl <- bal7
  bal7tl$tip.label <- bal7$tip[map]
  
  expect_equal(SortTree(Preorder(as.phylo.HPart(hbMap))), SortTree(bal7tl))
})
