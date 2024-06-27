library("TreeTools", quiet = TRUE)
if(!exists("pv")) pv <- function (x) x

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

test_that("confusion()", {
  skip(37)
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

test_that("SPR calculated correctly", {
  Tree <- function (txt) ape::read.tree(text = txt)
  expect_equal(.SPRPair(ape::read.tree(text = "((a, b), (c, d));"),
                        ape::read.tree(text = "((a, c), (b, d));")),
               1L)
  expect_equal(.SPRPair(PectinateTree(letters[1:26]),
                        PectinateTree(letters[c(2:26, 1)])),
               1L)
  expect_equal(.SPRPair(
    tree1 <- PectinateTree(letters[1:26]),
    tree2 <- Tree("(g, (h, (i, (j, (k, (l, ((m, (c, (b, a))), (n, (o, (p, (q, (r, (s, (t, (u, (v, (w, (x, (y, (z, (f, (e, d))))))))))))))))))))));")),
    2)
  expect_equal(.SPRPair(
    tree1 <- PectinateTree(letters[1:26]),
    tree2 <- Tree("(g, (h, (i, (j, (k, (l, (m, (n, (o, (p, (q, (r, (s, (t, (u, (v, (w, (x, (y, (z, (f, ((e, (c, (b, a))), d))))))))))))))))))))));")),
    2)
  
  expect_equal(.phangornSPRDist(PectinateTree(letters[1:26]),
                                PectinateTree(letters[c(2:26, 1)])),
               c(spr = 1L))
  
  
  set.seed(0)
  tr <- vector("list", 13)
  tr[[1]] <- Postorder(RandomTree(25, root = TRUE))
  for (i in seq_len(12) + 1L) {
    tr[[i]] <- Postorder(TreeSearch::SPR(tr[[i - 1]]))
  }
  
  expect_gt(.phangornSPRDist(tr[[3]], tr[[11]]), 7) # it's 8
  expect_equal(SPRDist(tr[[3]], tr[[11]]), 7)
  # expect_equal(TBRDist::USPRDist(tr[[3]], tr[[11]]), 7)
  
  if (interactive()) {
  #  trueDist <- TBRDist::USPRDist(tr)
    trueDist <- readRDS("true-25tip-12spr.Rds")
  
    
    
    par(mfrow = c(1, 2))
    distRange <- c(simDist - phanDist, simDist - bestDist)
    hist(distRange, col = NA, border = NA)
    hist(simDist - phanDist, add = TRUE, col = 2)
    hist(simDist - bestDist, add = TRUE, col = "#88ee4488")
    
    plot(simDist, simDist, type = "n", asp = 1, ylim = range(distRange),
         xlab = "Number of SPR moves")
    abline(0, 0, col = 3)
    jd <- jitter(simDist)
    #points(jd, trueDist, pch = 7, col = 3)
    #points(jd, phanDist, pch = 1)
    #points(jd, bestDist, pch = 3, col = 2)
    points(jd, phanDist - trueDist, pch = 5, col = 4)
    points(jd, bestDist - trueDist, pch = 4, col = 5)
  }
  
  nTip <- 130
  nSPR <- 35
  
  set.seed(0)
  tr <- vector("list", nSPR + 1L)
  tr[[1]] <- Postorder(RandomTree(nTip, root = TRUE))
  expect_equal(SPRDist(tr[[1]], tr[[1]]), 0)
  for (i in seq_len(nSPR) + 1L) {
    tr[[i]] <- Postorder(TreeSearch::SPR(tr[[i - 1]]))
  }
  
  phanDist <- .phangornSPRDist(tr)
  
  SPRDist(tr[[1]], tr)
  
  
  testDist <- SPRDist(tr)
  simDist <- dist(seq_along(tr))
  
  expect_true(all(testDist >= phanDist))
  expect_true(all(testDist <= simDist))
  
  # bestDist <- as.dist(pmin(as.matrix(testDist), as.matrix(SPRDist(rev(tr)))[rev(seq_along(tr)), rev(seq_along(tr))]))
  bestDist <- testDist # assert symmetry
  
  overShot <- as.matrix(testDist) > as.matrix(simDist)
  overs <- colSums(overShot) > 0
  overShot[overs, overs]
  
  underShot <- as.matrix(testDist) < as.matrix(phanDist)
  unders <- colSums(underShot) > 0
  underShot[unders, unders]
  
  tree1 <- tr[[1]]
  tree2 <- tr[[36]]
  .SPRPair(tree1, tree2, debug = TRUE)
  
  
  tree1 <- tr[[3]]
  tree2 <- tr[[11]]
  .SPRPair(tree1, tree2, debug = TRUE)
  
  tree1 <- tr[[14]]
  tree2 <- tr[[24]]
  .SPRPair(tree1, tree2, debug = TRUE)

  # ub(SPRDist(tr), .phangornSPRDist(tr), times = 3)
  # pv(testDist <- SPRDist(tr))

  
  if (interactive()) {
    skip("This shouldn't run!")
    if (nTip < 51 && nSPR < 13) {
      if (nTip == 25 && nSPR == 12) {
        trueDist <- readRDS("true-25tip-12spr.Rds")
      } else {
        trueDist <- TBRDist::USPRDist(tr)
      }
    }
  } else {
    trueDist <- simDist
  }
  
  
  par(mfrow = c(1, 2))
  distRange <- c(simDist - phanDist, simDist - bestDist)
  hist(distRange, col = NA, border = NA)
  hist(simDist - phanDist, add = TRUE, col = 2)
  hist(simDist - bestDist, add = TRUE, col = "#88ee4488")
  
  plot(simDist, simDist, type = "n", asp = 1, ylim = range(distRange),
       xlab = "Number of SPR moves")
  abline(0, 0, col = 3)
  jd <- jitter(simDist)
  #points(jd, trueDist, pch = 7, col = 3)
  #points(jd, phanDist, pch = 1)
  #points(jd, bestDist, pch = 3, col = 2)
  points(jd, phanDist - trueDist, pch = 5, col = 4)
  points(jd, bestDist - trueDist, pch = 4, col = 5)
})

test_that("SPR.dist called safely", {
  skip("#TODO")
  library("TreeTools")
  PhangornSPR <- function(...) unname(phangorn::SPR.dist(...))
  tree1 <- as.phylo(0, 6)
  tree2 <- BalancedTree(6)
  expect_equal(SPRDist(as.phylo(0, 6), BalancedTree(6)),
               PhangornSPR(Postorder(as.phylo(0, 6)),
                           Postorder(BalancedTree(6))))
  expect_equal(SPRDist(as.phylo(0:5, 6), BalancedTree(6)), 
               PhangornSPR(structure(lapply(as.phylo(0:5, 6), Postorder),
                                     class = "multiPhylo"),
                           Postorder(BalancedTree(6))))
  expect_equal(SPRDist(BalancedTree(6), as.phylo(0:5, 6)),
               SPRDist(as.phylo(0:5, 6), BalancedTree(6)))
  expect_equal(SPRDist(BalancedTree(6), PectinateTree(6)),
               SPRDist(list(BalancedTree(6), PectinateTree(6)))[[1]],
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
  trs12 <- structure(list(tr1, tr2), class = "multiPhylo")
  trs123 <- structure(list("one" = tr1, "two" = tr2, "thr" = tr3),
                      class = "multiPhylo")
  SprpS <- function(...) SPRDist(..., symmetric = TRUE)
  expect_equal(SprpS(tr1, tr3), SprpS(tr3, tr1))
  expect_equal(2L, length(SprpS(trs12, tr3)))
  expect_equal(SprpS(trs12, tr3), SprpS(tr3, trs12))
  expect_equal(SprpS(trs12, trs123), t(SprpS(trs123, trs12)))
  expect_equal(SprpS(rev(trs123), trs123), t(SprpS(trs123, rev(trs123))))
  expect_equal(1 - diag(1, 3), as.matrix(SprpS(trs123)),
               ignore_attr = TRUE)
})
