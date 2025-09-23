library("TreeTools")

test_that("as.phylo.HPart", {
  bal7 <- BalancedTree(7)
  hb7 <- as.HPart(bal7)
  expect_equal(Preorder(as.phylo.HPart(hb7)), bal7)
  
  bal17 <- BalancedTree(17)
  hb17 <- as.HPart(bal17)
  expect_equal(Preorder(as.phylo.HPart(hb17)), bal17)
})

test_that("as.HPart.numeric", {
  hpList <- as.HPart(list(list(1, 3, 9),
                          list(2, 4, 8),
                          list(5, 6, 7)))
  expect_equal(class(hpList), "HPart")
  
  hpNum <- as.HPart(c(1, 2, 1, 2, 3, 3, 3, 2, 1))
  expect_equal(class(hpNum), "HPart")
  
  expect_equal(SelfHMI(hpNum), SelfHMI(hpList))
  expect_equal(HMI(hpNum, hpList), SelfHMI(hpNum))
  
  flatP <- as.HPart(list(as.list(1:5), as.list(6:9)))
  hp9 <- as.HPart(BalancedTree(1:9))
  expect_equal(HMI(flatP, hp9), 0.99107606)
  
})

test_that("as.HPart.unimplemented", {
  expect_error(as.HPart(matrix()), "no applicable method")
  expect_error(as.HPart(list(letters, LETTERS)), "leaves must be integers")
  expect_error(as.HPart(list(list(1, 2, 3), list(0, 1, 2))),
               ".eaves must contain each integer")
  expect_error(as.HPart(list(list(1, 2, 3), list(3, 1, 2))),
               ".eaves must contain each integer")
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

test_that("plot.HPart", {
  skip_if_not_installed("vdiffr")
  vdiffr::expect_doppelganger("plot-HPart", function() 
    plot(as.HPart(list(list(1, 2, 3), list(4, list(5, 6)))))
  )
})

test_that("Renumber.HPart", {
  expect_error(RenumberTips(as.HPart(list(1, 2, 4, 3)), 4:2),
               "labels 1 missing from `tipOrder`")
})
