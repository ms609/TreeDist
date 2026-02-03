test_that("SPR: keep_and_reroot()", {
  library("TreeTools", quietly = TRUE)
  
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
  library("TreeTools", quietly = TRUE)
  
  expect_error(
    mismatch_size(as.Splits(c(T, T, F)), as.Splits(c(T, T, T, T))),
    "differ in `nTip"
  )
  expect_error(
    mismatch_size(matrix(as.raw(3), 1, 1), as.Splits(c(T, T, T, T))),
    "nTip attribute"
  )
  expect_error(
    mismatch_size(as.Splits(c(T, T, T, T)), matrix(as.raw(3), 1, 1)),
    "nTip attribute"
  )
  expect_error(
    mismatch_size(as.Splits(matrix(T, 2, 4)), as.Splits(c(T, T, T, T))),
    "number of splits"
  )
  splits <- as.Splits(rbind(c(T, T, T, F, F), c(T, F, F, F, T)))
  Test <- function(s1, s2) {
    expect_equal(length(s1), length(s2))
    nSplits <- length(s1)
    i <- rep(seq_len(nSplits), nSplits)
    j <- rep(seq_len(nSplits), each = nSplits)
    expect_equal(
      mismatch_size(s1, s2),
      TipsInSplits(xor(s1[[i]], s2[[j]]), smallest = TRUE)
    )
  }
  Test(as.Splits(c(T, T, T, F, F)), as.Splits(c(T, F, F, F, T)))

  set.seed(0)
  splits <- as.Splits(t(replicate(10, sample(c(T, F), 99, replace = TRUE))))
  Test(splits[[1]], splits[[2]])
  Test(splits[[1:2]], splits[[2:3]])
  Test(splits, rev(splits))
})

test_that("confusion()", {
  library("TreeTools", quietly = TRUE)
  
  TestConfusion <- function(a, b) {
    i <- rep(seq_along(a), each = length(b))
    j <- rep(seq_along(b), length(a))
    expect_equal(
      confusion(a, b),
      aperm(array(
        c(TipsInSplits(a[[i]] & b[[j]]),
          TipsInSplits(a[[i]] & !b[[j]]),
          TipsInSplits(!a[[i]] & b[[j]]),
          TipsInSplits(!a[[i]] & !b[[j]])
        ),
        c(length(a), length(b), 4)),
        c(3, 2, 1)
      )
    )
  }

  TestConfusion(as.Splits(c(T, T, T, F, F)), as.Splits(c(T, F, F, F, T)))

  set.seed(0)
  splits <- as.Splits(t(replicate(10, sample(c(T, F), 99, replace = TRUE))))
  TestConfusion(splits[[1]], splits[[2]])
  TestConfusion(splits[[1:2]], splits[[2:3]])
  TestConfusion(splits, rev(splits))
})

test_that("SPR deOliveira2008 calculation looks valid", {
  library("TreeTools", quietly = TRUE)
  
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
  tr <- vector("list", nSPR + 1L)
  tr[[1]] <- Postorder(TreeTools::RandomTree(nTip, root = TRUE))
  expect_equal(SPRDist(tr[[1]], tr[[1]])[[1]], 0)
  for (i in seq_len(nSPR) + 1L) {
    tr[[i]] <- Postorder(TreeSearch::SPR(tr[[i - 1]]))
  }
  
  testDist <- as.matrix(SPRDist(tr, method = "de Oliv"))
  simDist <- as.matrix(dist(seq_along(tr)))
  for (i in 1:nSPR) for (j in 2:nSPR) {
    {if (i < j) expect_gte else expect_lte}(testDist[i, j], testDist[i, j - 1])
  }
  
  overShot <- as.matrix(testDist) > as.matrix(simDist)
  overs <- colSums(overShot) > 0
  overShot[overs, overs]
  expect_false(any(overs))
})

test_that("SPR calculated correctly", {
  library("TreeTools", quietly = TRUE)
  
  expect_exact <- function(x, y, method = "confl") {
    if (is.character(x)) x <- Tree(x)
    if (is.character(y)) y <- Tree(y)
    expect_equal(
      SPRDist(x, y, method = method)[[1]],
      TBRDist::USPRDist(x, y)
    )
  }
  
  Tree <- function(txt) ape::read.tree(text = txt)

  expect_equal(
    .SPRConfl(
      ape::read.tree(text = "((a, b), (c, d));"),
      ape::read.tree(text = "((a, c), (b, d));")
    )[[1]],
    1L
  )
  expect_equal(
    .SPRConfl(PectinateTree(letters[1:26]), PectinateTree(letters[c(2:26, 1)]))[[1]],
    1L
  )
  
  options(sprH = "viNorm")
  # Looks simple, but requires ties to be broken suitably
  # Passes with conf, viNorm (needed tiebreaker in each case)
  # Fails with joint, ami, vi: tiebreaker not yet implemented!
  expect_exact("(a,(d,(b,(c,X))));", "(a,((b,c),(X,d)));") # distance = 1
  
  expect_exact("((((b,c),d),e),a);", "(a,(b,((e,c),d)));") # distance = 2
  expect_failure(
    expect_exact("((((b,c),d),e),a);", "(a,(b,((e,c),d)));", "deo")
  )
  
  # Example AZ33: IJK and DEF are schlepped
  # Passes with joint, ami, viNorm
  # Fails with vi (5), conf (7)
  expect_equal(
    .SPRConfl(
      tree1 <- PectinateTree(letters[1:26]),
      tree2 <- Tree(
        "(g, (h, (i, (j, (k, (l, ((m, (c, (b, a))), (n, (o, (p, (q, (r, (s, (t, (u, (v, (w, (x, (y, (z, (f, (e, d))))))))))))))))))))));"
      )
    )[[1]],
    2
  )
  
  # Benefits from the "divide and conquer" step, 2026-02-03
  expect_equal(
    .SPRConfl(
      tree1 <- PectinateTree(letters[1:26]),
      tree2 <- Tree(
        "(g, (h, (i, (j, (k, (l, (m, (n, (o, (p, (q, (r, (s, (t, (u, (v, (w, (x, (y, (z, (f, ((e, (c, (b, a))), d))))))))))))))))))))));"
      )
    )[[1]],
    2
  )
  
  lockedMid1 <- Tree("((((a1, a2), a3), ((b1, b2), b3)),
                     (((c1, c2), c3), ((d1, d2), d3)));")
  lockedMid2 <- Tree("(((a1, (a2, a3)), (c1, (c2, c3))),
                     ((b1, (b2, b3)), (d1, (d2, d3))));")
  expect_equal(.SPRConfl(lockedMid1, lockedMid2)[[1]], 5)

  set.seed(0)
  tr <- vector("list", 13)
  tr[[1]] <- Postorder(RandomTree(25, root = TRUE))
  for (i in seq_len(12) + 1L) {
    tr[[i]] <- Postorder(TreeSearch::SPR(tr[[i - 1]]))
  }

  write.tree(tr[[3]])
  
  expect_equal(SPRDist(tr[[3]], tr[[11]], method = "DeO"), 8)
  # 11 with divide and conquer (9 with tie-breaker); 9 without
  expect_equal(SPRDist(tr[[3]], tr[[11]], method = "confl"), 7) #
  expect_equal(SPRDist(tr[[3]], tr[[11]], method = "exp"), 7) #
  # expect_equal(TBRDist::USPRDist(tr[[3]], tr[[11]]), 7) # at 2026-01-15: Actually 8
  
  
  # Simplified example for reproducibility
  Simplify <- function(tr) {
    # Critical to the behaviour: t19, t25, t5
    tr |>
      DropTip("t17") |> # Reduces by a step
      DropTip("t13") |> # Reduces by a step
      DropTip("t8") |> # Reduces by a step
      DropTip("t2") |> # Reduces by a step
      DropTip("t19") |> # Reduces by a step
      DropTip(c("t1", "t4", "t6", "t11", "t12", "t14", "t15", "t9", "t24")) # No difference to score
  }
  # t3 <- Simplify(tr[[3]])
  # t11 <- Simplify(tr[[11]])
  t3 <- Tree("((((t21,(((((t5,t22),t7),t20),t25),(t23,t16))),t18),t10),t3);")
  t11 <- Tree("((((t21,t18),t10),(((t20,t5),(t7,t22)),(t23,(t16,t25)))),t3);")
  
  deO <- SPRDist(t3, t11, method = "DeO")
  conf <- SPRDist(t3, t11, method = "confl")
  
  expect_equal(conf - deO,
               SPRDist(tr[[3]], tr[[11]], method = "confl") -
                 SPRDist(tr[[3]], tr[[11]], method = "DeO"))
  
  

  nTip <- 130
  nSPR <- 35

  set.seed(0)
  tr <- vector("list", nSPR + 1L)
  tr[[1]] <- Postorder(RandomTree(nTip, root = TRUE))
  expect_equal(SPRDist(tr[[1]], tr[[1]])[[1]], 0)
  for (i in seq_len(nSPR) + 1L) {
    tr[[i]] <- Postorder(TreeSearch::SPR(tr[[i - 1]]))
  }

  phanDist <- SPRDist(tr, method = "deO")

  SPRDist(tr[[1]], tr)

  testDist <- SPRDist(tr)
  simDist <- dist(seq_along(tr))

  expect_true(all(testDist <= simDist))

  # bestDist <- as.dist(pmin(as.matrix(testDist), as.matrix(SPRDist(rev(tr)))[rev(seq_along(tr)), rev(seq_along(tr))]))
  bestDist <- testDist # assert symmetry

  overShot <- as.matrix(testDist) > as.matrix(simDist)
  overs <- colSums(overShot) > 0
  overShot[overs, overs]

  underShot <- as.matrix(testDist) < as.matrix(phanDist)
  unders <- colSums(underShot) > 0
  underShot[unders, unders]

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

  expect_true(all(testDist >= phanDist))

  # Cases to debug
  opt <- options(debugSPR = TRUE)
  on.exit(options(opt))
  tree1 <- tr[[1]]
  tree2 <- tr[[36]]
  .SPRConfl(tree1, tree2)

  tree1 <- tr[[3]]
  tree2 <- tr[[11]]
  .SPRConfl(tree1, tree2)

  tree1 <- tr[[14]]
  tree2 <- tr[[24]]
  .SPRConfl(tree1, tree2)
  options(opt)

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
