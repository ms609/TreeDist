library(ape)
library(phangorn)

context("tree scoring")
comb11 <- read.tree(text="(a, (b, (c, (d, (e, (f, (g, (h, (i, (j, k))))))))));")
data11 <- cbind(upper.tri(matrix(FALSE, 11, 11))[, 3:10], lower.tri(matrix(FALSE, 11, 11))[, 2:9])
rownames(data11) <- letters[1:11]
phy11 <- phyDat(data11, type='USER', levels=c(FALSE, TRUE))

test_that("Tree is scored correctly", {
  expect_equal(FitchScore(comb11, phy11), 16)
  expect_equal(FitchScore(Pruningwise(comb11), phy11), 16)
  expect_equal(FitchScore(Cladewise(comb11), phy11), 16)
  
  expect_equal(FitchScore(comb11, data11, TipData=TipsAreRows), 16)
})


### # TODO : TipsAreRows doesn't send the data in the right way.
###
###    FitchFunction = C_Fitch_Score
###  treeOrder <- attr(tree, 'order')
###  if (is.null(treeOrder) || treeOrder != 'postorder') tree <- Postorder(tree)
###  treeEdge <- tree$edge
###  parent <- treeEdge[, 1]
###  child <- treeEdge[, 2]
###  tipLabel <- tree$tip.label
###  nr <- at$nr
###  if (is.null(at$nr)) {
###    nr <- dim(data)
###    nr <- if (nr[1] == length(tipLabel)) nr[2] else nr[1]
###  }
###  charWeights <- at$weight
###  if (is.null(charWeights)) charWeights <- rep(1, nr)
###  if (class(data) =='phyDat') {
###    levs <- attr(data, 'levels')
###    contrast <- attr(data, 'contrast')
###    index <- as.integer(contrast %*% 2L ^ (seq_along(attr(data, 'levels')) - 1))
###    reformedData <- vapply(data, function (X) index[X], integer(nr))
###    characters <- TipsAreColumns(reformedData, tipLabel)
###  } else {
###    characters <- TipData(data, tipLabel)
###  }
###  
###  C_Fitch(
###      characters = characters, 
###      nChar = nr,
###      parent, child,
###      nEdge = length(parent), 
###      weight = charWeights, 
###      maxNode = max(parent), #parent[1] IF tree in postorder
###      nTip = length(tipLabel)
###    )
###  