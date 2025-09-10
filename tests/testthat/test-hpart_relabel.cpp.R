library("TreeTools")

test_that("as.phylo.HPart", {
  bal7 <- BalancedTree(7)
  plot(bal7)
  hb7 <- as.HPart_cpp(bal7)
  expect_equal(Preorder(as.phylo.HPart_cpp(hb7)), bal7)
  
  bal17 <- BalancedTree(17)
  plot(bal17)
  hb17 <- as.HPart_cpp(bal17)
  expect_equal(Preorder(as.phylo.HPart_cpp(hb17)), bal17)
})

test_that("HParts are relabelled correctly", {
  bal7 <- BalancedTree(7)
  hb7 <- as.HPart_cpp(bal7)
  
  map <- c(7:4, 1:3)
  mappedLabels <- paste0("t", map)
  
  hbMap <- RenumberTips(hb7, mappedLabels)
  # Here we want only to map the internal node IDs
  attr(hbMap, "tip.label") <- TipLabels(hb7)
  
  bal7tl <- bal7
  bal7tl$tip.label <- bal7$tip[map]
  
  expect_equal(SortTree(Preorder(as.phylo.HPart_cpp(hbMap))), SortTree(bal7tl))
})
