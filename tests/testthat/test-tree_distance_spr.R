library("TreeTools", quiet = TRUE, warn.conflicts = FALSE)
if(!exists("pv")) pv <- function (x) x

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
    expect_equal(mismatch_size(s1, s2), TipsInSplits(xor(s1[[i]], s2[[j]])))
  }
  Test(as.Splits(c(T, T, T, F, F)), as.Splits(c(T, F, F, F, T)))
  
  set.seed(0)
  splits <- as.Splits(t(replicate(10, sample(c(T, F), 99, replace = TRUE))))
  Test(splits[[1]], splits[[2]])
  Test(splits[[1:2]], splits[[2:3]])
  Test(splits, rev(splits))
})

test_that("SPR calculated correctly", {
  expect_equal(.SPRPair(ape::read.tree(text = "((a, b), (c, d));"),
                        ape::read.tree(text = "((a, c), (b, d));")),
               1L)
  
  nTip <- 90
  nTrees <- 1
  nSPR <- 30

  set.seed(0)
  tr <- vector("list", nSPR + 1L)
  tr[[1]] <- Postorder(RandomTree(nTip, root = TRUE))
  for (i in seq_len(nSPR) + 1L) {
    tr[[i]] <- Postorder(TreeSearch::SPR(tr[[i - 1]]))
  }
  pv(testDist <- SPRDist(tr))
  bestDist <- as.dist(pmax(as.matrix(testDist), as.matrix(SPRDist(rev(tr)))[rev(seq_len(nSPR + 1)), rev(seq_len(nSPR + 1))]))
  expect_equal(SPRDist(tr[[1]], tr[[1]]), 0)
  
  trueDist <- dist(seq_len(nSPR + 1) - 1)
  # ub(SPRDist(tr), .phangornSPRDist(tr), times = 4)
  phanDist <- .phangornSPRDist(tr)

  par(mfrow = c(1, 2))  
  hist(trueDist - phanDist, col = 2)
  hist(trueDist - bestDist, add = TRUE, col = "#88ee4488")
  
  plot(trueDist, trueDist, type = "n", asp = 1,
       xlab = "Number of SPR moves")
  abline(0, 1)
  jd <- jitter(trueDist)
  points(jd, phanDist, pch = 1)
  points(jd, bestDist, pch = 3, col = 2)
  arrows(jd, phanDist, jd, bestDist, col = 3)
  summary(lm(trueDist~phanDist + 0))
  summary(lm(trueDist~bestDist + 0))
})

test_that("SPR.dist called safely", {
  library("TreeTools")
  PhangornSPR <- function(...) unname(phangorn::SPR.dist(...))
  tree1<-as.phylo(0, 6)
  tree2=BalancedTree(6)
  expect_equal(SPRDist(as.phylo(0, 6), BalancedTree(6)),
               PhangornSPR(Postorder(as.phylo(0, 6)),
                           Postorder(BalancedTree(6))))
  expect_equal(SPRDist(as.phylo(0:5, 6), BalancedTree(6)), 
               PhangornSPR(structure(lapply(as.phylo(0:5, 6), Postorder),
                                     class = 'multiPhylo'),
                           Postorder(BalancedTree(6))))
  expect_equal(SPRDist(BalancedTree(6), as.phylo(0:5, 6)),
               SPRDist(as.phylo(0:5, 6), BalancedTree(6)))
  expect_equal(SPRDist(BalancedTree(6), PectinateTree(6)),
               SPRDist(list(BalancedTree(6), PectinateTree(6)))[1],
               ignore_attr = TRUE)
  
  # https://github.com/KlausVigo/phangorn/issues/97
  tr1 <- structure(list(edge = structure(c(11L, 11L, 10L, 10L, 9L, 9L, 8L, 8L, 7L,
                                           7L, 2L, 6L, 5L, 11L, 4L, 10L, 3L, 9L, 
                                           1L, 8L), .Dim = c(10L, 2L)),
                        tip.label = c("t1", "t2", "t3", "t4", "t5", "t6"), 
                        Nnode = 5),
                   class = "phylo")
  
  tr2 <- structure(list(edge = structure(c(10L, 10L, 11L, 11L, 9L, 9L, 8L, 8L, 7L,
                                           7L, 2L, 6L, 3L, 4L, 5L, 10L, 11L, 9L, 
                                           1L, 8L), .Dim = c(10L, 2L)),
                        tip.label = c("t1", "t2", "t3", "t4", "t5", "t6"),
                        Nnode = 5),
                   class = "phylo")
  
  tr3 <- structure(list(edge = structure(c(9L, 9L, 11L, 11L, 10L, 10L, 8L, 8L, 7L, 
                                           7L, 1L, 2L, 4L, 5L, 6L, 11L, 3L, 9L, 
                                           8L, 10L), .Dim = c(10L, 2L)),
                        Nnode = 5L, 
                        tip.label = c("t1", "t2", "t3", "t4", "t5", "t6")),
                   class = "phylo")
  trs12 <- structure(list(tr1, tr2), class = 'multiPhylo')
  trs123 <- structure(list('one' = tr1, 'two' = tr2, 'thr' = tr3),
                      class = 'multiPhylo')
  SprpS <- function(...) SPRDist(..., symmetric = TRUE)
  expect_equal(SprpS(tr1, tr3), SprpS(tr3, tr1))
  expect_equal(2L, length(SprpS(trs12, tr3)))
  expect_equal(SprpS(trs12, tr3), SprpS(tr3, trs12))
  expect_equal(SprpS(trs12, trs123), t(SprpS(trs123, trs12)))
  expect_equal(SprpS(rev(trs123), trs123), t(SprpS(trs123, rev(trs123))))
  expect_equal(1 - diag(1, 3), as.matrix(SprpS(trs123)),
               ignore_attr = TRUE)
})
