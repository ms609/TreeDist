library("TreeTools")

test_that("HParts are relabelled correctly", {
  bal7 <- BalancedTree(7)
  plot(bal7)
  hb7 <- as.HPart_cpp(bal7)
  expect_equal(Preorder(as.phylo.HPart_cpp(hb7)), bal7)
  
  bal17 <- BalancedTree(17)
  plot(bal17)
  hb17 <- as.HPart_cpp(bal17)
  expect_equal(Preorder(as.phylo.HPart_cpp(hb17)), bal17)
})
