library("TreeTools", quietly = TRUE)

test_that("SPR: keep_and_reroot()", {
  tree1 <- Postorder(BalancedTree(12))
  tree2 <- Postorder(PectinateTree(12))
  keep <- as.logical(tabulate(8:12, 12))

  result <- keep_and_reroot(tree1, tree2, keep)
  expect_equal(result[[1]], RootTree(KeepTip(tree1, keep), 1))
  expect_equal(result[[2]], RootTree(KeepTip(tree2, keep), 1))
  
  reduced <- keep_and_reduce(tree1, tree2, keep)
  expect_equal(Preorder(reduced[[1]]), Preorder(DropTip(result[[1]], "t9")))
  expect_equal(Preorder(reduced[[2]]), Preorder(DropTip(result[[2]], "t9")))
  
  skip_if_not_installed("TreeTools", "1.11.1.9002")
  twoZeroes <- list(Preorder(ZeroTaxonTree()), Preorder(ZeroTaxonTree()))
  expect_equal(keep_and_reroot(SingleTaxonTree(), SingleTaxonTree(), FALSE),
               twoZeroes)
  expect_equal(keep_and_reduce(SingleTaxonTree(), SingleTaxonTree(), TRUE),
               twoZeroes)
})

test_that("SPR: Under the hood", {
  expect_error(mismatch_size(as.Splits(c(T, T, F)), as.Splits(c(T, T, T, T))),
               "differ in `nTip")
  expect_error(mismatch_size(matrix(as.raw(3), 1, 1), 
                             as.Splits(c(T, T, T, T))),
               "nTip attribute")
  expect_error(mismatch_size(as.Splits(c(T, T, T, T)),
                             matrix(as.raw(3), 1, 1)),
               "nTip attribute")
  expect_error(mismatch_size(as.Splits(matrix(T, 2, 4)),
                             as.Splits(c(T, T, T, T))),
               "number of splits")
  splits <- as.Splits(rbind(c(T, T, T, F, F),
                            c(T, F, F, F, T)))
  Test <- function (s1, s2) {
    expect_equal(length(s1), length(s2))
    nSplits <- length(s1)
    i <- rep(seq_len(nSplits), nSplits)
    j <- rep(seq_len(nSplits), each = nSplits)
    expect_equal(mismatch_size(s1, s2),
                 TipsInSplits(xor(s1[[i]], s2[[j]]), smallest = TRUE))
  }
  Test(as.Splits(c(T, T, T, F, F)), as.Splits(c(T, F, F, F, T)))
  
  set.seed(0)
  splits <- as.Splits(t(replicate(10, sample(c(T, F), 99, replace = TRUE))))
  Test(splits[[1]], splits[[2]])
  Test(splits[[1:2]], splits[[2:3]])
  Test(splits, rev(splits))
})

test_that("confusion() fails gracefully", {
  x <- as.Splits(c(T, T, T, F, F))
  xNoTip <- x
  attr(xNoTip, "nTip") <- NULL
  xx <- as.Splits(rep(c(T, F, F, F, T), 4))
  expect_error(confusion(x, xx), "differ in `nTip`")
  expect_error(confusion(xNoTip, x), "`x` lacks nTip attribute")
  expect_error(confusion(x, xNoTip), "`y` lacks nTip attribute")
  expect_error(confusion(x, xx), "differ in `nTip`")
  expect_error(confusion(c(x, x), x), "number of splits")
})

test_that("confusion()", {
  TestConfusion <- function (a, b) {
    i <- rep(seq_along(a), each = length(b))
    j <- rep(seq_along(b), length(a))
    expect_equal(
      confusion(a, b),
      aperm(array(c(TipsInSplits(a[[i]] & b[[j]]),
                    TipsInSplits(a[[i]] & !b[[j]]),
                    TipsInSplits(!a[[i]] & b[[j]]),
                    TipsInSplits(!a[[i]] & !b[[j]])),
                  c(length(a), length(b), 4)), c(3, 2, 1))
    )
  }
  
  TestConfusion(as.Splits(c(T, T, T, F, F)), as.Splits(c(T, F, F, F, T)))
  
  set.seed(0)
  splits <- as.Splits(t(replicate(10, sample(c(T, F), 99, replace = TRUE))))
  TestConfusion(splits[[1]], splits[[2]])
  TestConfusion(splits[[1:2]], splits[[2:3]])
  TestConfusion(splits, rev(splits))
})

test_that("SPRDist handles input formats", {
  bal9 <- BalancedTree(9)
  pec9 <- PectinateTree(9)
  
  expect_null(SPRDist(bal9))
  expect_equal(SPRDist(bal9, bal9), 0)
  
  expect_equal(SPRDist(list(bal9, bal9), bal9), c(0, 0))
  expect_equal(SPRDist(c(bal9, bal9), bal9), c(0, 0))
  expect_equal(SPRDist(bal9, list(bal9, bal9)), c(0, 0))
  expect_equal(SPRDist(bal9, c(bal9, bal9)), c(0, 0))
  
  expect_equal(SPRDist(list(bal9, pec9, pec9), list(pec9, bal9)),
               matrix(c(2, 0, 0, 0, 2, 2), 3, 2))
  expect_equal(SPRDist(c(bal9, pec9, pec9), c(pec9, bal9)),
               matrix(c(2, 0, 0, 0, 2, 2), 3, 2))
  
  self <- SPRDist(list(bal9, pec9))
  at <- attributes(self)
  expect_equal(at[["Size"]], 2)
  expect_equal(at[["class"]], "dist")
  expect_equal(at[["Diag"]], FALSE)
  expect_equal(at[["Upper"]], TRUE)
  expect_equal(self[[1]], dist(c(0, 2))[[1]])
})

test_that("SPR deOliveira2008 calculation looks valid", {
  # We do not expect to obtain identical results to phangorn::SPR.dist,
  # because ties are broken in a different arbitrary manner.
  # We're thus left with quite a loose test.
  Tree <- function (txt) ape::read.tree(text = txt)
  
  expect_equal(SPRDist(PectinateTree(letters[1:26]),
                       PectinateTree(letters[c(2:26, 1)]),
                       method = "deOliv"),
               1L)
  
  nTip <- 130
  nSPR <- 35
  
  set.seed(0)
  skip_if_not_installed("TreeSearch")
  tr <- vector("list", nSPR + 1L)
  tr[[1]] <- Postorder(TreeTools::RandomTree(nTip, root = TRUE))
  expect_equal(SPRDist(tr[[1]], tr[[1]]), 0)
  for (i in seq_len(nSPR) + 1L) {
    tr[[i]] <- Postorder(TreeSearch::SPR(tr[[i - 1]]))
  }
  
  testDist <- as.matrix(SPRDist(tr, method = "de Oliv"))
  simDist <- as.matrix(dist(seq_along(tr)))
  biggerThyNeighbour <- sapply(1:nSPR, function(i) sapply(2:nSPR, function(j)
    if (i < j) {
      testDist[i, j] - testDist[i, j - 1]
    } else {
      testDist[i, j - 1] - testDist[i, j]
    }
  ))
  
  errors <- biggerThyNeighbour < 0
  # We may see a few "errors" due to chance, but expect these to be rare
  rare <- nSPR
  expect_lt(sum(errors), rare)
  if (interactive() && any(errors)) {
    testDist[colSums(errors) > 0, ]
  }
  
  overShot <- as.matrix(testDist) > as.matrix(simDist)
  # We may overshoot where there are "knots" in trees and the optimal SPR
  # path is not equivalent to just pruning shared subtrees; these cases
  # ought to be rare.
  rare <- 0.10
  expect_lt(sum(overShot) / length(overShot), rare)
  
  if (interactive()) {
    # View these cases:
    overs <- colSums(overShot) > 0
    overShot[overs, overs]
  }
})
