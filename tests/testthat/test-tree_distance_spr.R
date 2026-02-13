test_that("SPR: keep_and_reroot()", {
  library("TreeTools", quietly = TRUE)
  
  tree1 <- Postorder(BalancedTree(12))
  tree2 <- Postorder(PectinateTree(12))
  keep <- as.logical(tabulate(8:12, 12))
  
  result <- keep_and_reroot(tree1, tree2, keep)
  expect_equal(result[[1]], RootTree(KeepTip(tree1, keep), 1))
  expect_equal(result[[2]], RootTree(KeepTip(tree2, keep), 1))
  
  result1 <- keep_and_reroot(tree1, tree2, 1:12 %in% 4)
  expect_equal(result1[[1]], Preorder(SingleTaxonTree("t4")))
  expect_equal(result1[[2]], Preorder(SingleTaxonTree("t4")))
  
  result0 <- keep_and_reduce(tree1, tree2, logical(12))
  expect_equal(result0[[1]], Preorder(ZeroTaxonTree()))
  
  reduced <- keep_and_reduce(tree1, tree2, keep)
  expect_equal(Preorder(reduced[[1]]), Preorder(DropTip(result[[1]], "t9")))
  expect_equal(Preorder(reduced[[2]]), Preorder(DropTip(result[[2]], "t9")))
  
  result0 <- keep_and_reduce(tree1, tree2, 1:12 %in% 4)
  expect_null(result0[[1]])
})

test_that("SPRDist() handles different inputs", {
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

test_that("SPR shortcuts okay - shape-pairs", {
  library("TreeTools", quietly = TRUE)
  
  depth <- if (getOption("slowMode", FALSE)) 5 else 1
  test_spr_exact <- function(nTip, nRep = 1, seed = 0) {
    set.seed(seed)
    combs <- as.integer(NUnrootedShapes(nTip))
    nIter <- nRep * depth
    for (i in seq_len(combs) - 1) {
      t1 <- UnrootedTreeWithShape(i, nTip, tipLabels = TipLabels(nTip))
      for (j in seq_len(combs) - 1) {
        t2 <- UnrootedTreeWithShape(j, nTip)
        res <- vapply(seq_len(nIter), function(k) {
          t2$tip.label <- paste0("t", sample.int(nTip))
          c(.SPRRogue(t1, t2), TBRDist::USPRDist(t1, t2))
        }, double(2))
        expect_equal(res[1, ], res[2, ])
      }
    }
  }
  test_spr_exact(6, 20)
  test_spr_exact(7, 4)
})

test_that("SPR shortcuts okay - exhaustive", {
  skip_if_not(getOption("slowMode", FALSE))
  library("TreeTools", quietly = TRUE)
  nTip <- 7
  allTrees <- as.phylo(seq_len(NUnrooted(nTip)), nTip)
  
  sum(apply(combn(length(allTrees), 2), 2, function(ij) {
    reduced <- ReduceTrees(allTrees[[ij[[1]]]], allTrees[[ij[[2]]]])
    r1 <- reduced[[1]]
    if (is.null(r1) || NTip(r1) != nTip) return(NA)
    r2 <- reduced[[2]]
    exact <- TBRDist::USPRDist(r1, r2)
    shortcut <- SPRDist(r1, r2, method = "Rogue")
    equal <- exact == shortcut
    if (!equal) message("Mismatch: ", paste(ij, collapse = ", "))
    expect_true(equal)
    equal
  }), na.rm = TRUE)
})

test_that("SPR shortcuts okay - known answer", {
  library("TreeTools", quietly = TRUE)
  set.seed(0)
  trees <- lapply(1:8, function(XX) RandomTree(9, root = TRUE))
    
  exact <- apply(combn(length(trees), 2), 2, function(ij) {
    TBRDist::USPRDist(trees[[ij[[1]]]], trees[[ij[[2]]]])
  })
  opt <- options("sprShortcut" = 0)
  on.exit(options(opt))
  noCuts <- SPRDist(trees, method = "rogue")
  
  options("sprShortcut" = 4)
  cuts4 <- SPRDist(trees, method = "rogue")
  expect_true(all(cuts4 == noCuts))
  
  options("sprShortcut" = 5)
  cuts5 <- SPRDist(trees, method = "rogue")
  expect_true(all(cuts5 <= noCuts))
  expect_true(all(cuts5 >= exact))
  
  options("sprShortcut" = 6)
  cuts6 <- SPRDist(trees, method = "rogue")
  expect_true(all(cuts6 <= noCuts))
  expect_true(all(cuts6 >= exact))
  # Aspirational:
  expect_true(all(cuts6 == exact))
  
  options("sprShortcut" = 7)
  cuts7 <- SPRDist(trees, method = "rogue")
  expect_true(all(cuts7 <= noCuts))
  expect_true(all(cuts7 >= exact))
  # Aspirational:
  expect_true(all(cuts7 == exact))
})

test_that("SPR shortcuts okay - larger trees", {
  library("TreeTools", quietly = TRUE)
  set.seed(0)
  trees <- lapply(1:8, function(XX) RandomTree(45, root = TRUE))
  opt <- options("sprShortcut" = 0)
  on.exit(options(opt))
  noCuts <- SPRDist(trees, method = "rogue")
  
  options("sprShortcut" = 4)
  cuts4 <- SPRDist(trees, method = "rogue")
  expect_true(all(cuts4 == noCuts))
  
  options("sprShortcut" = 5)
  cuts5 <- SPRDist(trees, method = "rogue")
  expect_true(all(cuts5 <= noCuts))
  
  options("sprShortcut" = 6)
  cuts6 <- SPRDist(trees, method = "rogue")
  expect_true(all(cuts6 <= noCuts))
})

test_that("SPR calculated correctly", {
  library("TreeTools", quietly = TRUE)
  
  expect_exact <- function(x, y, method = "rogue") {
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
    .SPRRogue(
      ape::read.tree(text = "((a, b), (c, d));"),
      ape::read.tree(text = "((a, c), (b, d));")
    )[[1]],
    1L
  )
  expect_equal(
    .SPRConfl(PectinateTree(letters[1:26]), PectinateTree(letters[c(2:26, 1)]))[[1]],
    1L
  )
  expect_equal(
    .SPRRogue(PectinateTree(letters[1:26]), PectinateTree(letters[c(2:26, 1)]))[[1]],
    1L
  )
  
  options(sprH = "viNorm")
  options(sprH = "ami")
  # Looks simple, but requires ties to be broken suitably
  expect_exact("(a,(d,(b,(c,X))));", "(a,((b,c),(X,d)));", "confl") # distance = 1
  expect_exact("(a,(d,(b,(c,X))));", "(a,((b,c),(X,d)));", "rogue") # distance = 1
  
  expect_exact("((((b,c),d),e),a);", "(a,(b,((e,c),d)));", "confl") # distance = 2
  expect_exact("((((b,c),d),e),a);", "(a,(b,((e,c),d)));", "rogue")
  expect_failure(
    expect_exact("((((b,c),d),e),a);", "(a,(b,((e,c),d)));", "deo")
  )
  
  # Passes with ami, joint, vi
  # Fails with viNorm - should tiebreaker be non-normalized?
  expect_exact("(((t21,(((t15,t12),t23),t24)),t18),t19);",
               "((t21,t18),(((t12,(t19,t15)),t23),t24));", "confl")
  expect_exact("(((t21,(((t15,t12),t23),t24)),t18),t19);",
               "((t21,t18),(((t12,(t19,t15)),t23),t24));", "rogue")
  
  # Example AZ33: IJK and DEF are schlepped
  # Passes with joint, ami, viNorm
  # Fails with vi (5, 6 with tiebreaker), conf (7)
  expect_equal(
    .SPRConfl(
      tree1 <- PectinateTree(letters[1:26]),
      tree2 <- Tree(
        "(g, (h, (i, (j, (k, (l, ((m, (c, (b, a))), (n, (o, (p, (q, (r, (s, (t, (u, (v, (w, (x, (y, (z, (f, (e, d))))))))))))))))))))));"
      )
    )[[1]],
    2
  )
  
  expect_equal(
    .SPRRogue(
      tree1 <- PectinateTree(letters[1:26]),
      tree2 <- Tree(
        "(g, (h, (i, (j, (k, (l, ((m, (c, (b, a))), (n, (o, (p, (q, (r, (s, (t, (u, (v, (w, (x, (y, (z, (f, (e, d))))))))))))))))))))));"
      )
    )[[1]],
    2
  )
  
  # Requires "divide and conquer" step
  expect_equal(
    .SPRConfl(
      tree1 <- PectinateTree(letters[1:26]),
      tree2 <- Tree(
        "(g, (h, (i, (j, (k, (l, (m, (n, (o, (p, (q, (r, (s, (t, (u, (v, (w, (x, (y, (z, (f, ((e, (c, (b, a))), d))))))))))))))))))))));"
      )
    )[[1]],
    2
  )
  
  expect_equal(
    .SPRRogue(
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
  expect_equal(.SPRRogue(lockedMid1, lockedMid2)[[1]], 5)
  
  # Aspirational
  expect_exact("((((b3,b2),b1),(((d2,d1),((e3,e2),e1)),c)),a);",
               "((((d2,e3),d1),c),(((((e2,b3),e1),b1),b2),a));", method = "Rogue")

  set.seed(0)
  tr <- vector("list", 13)
  tr[[1]] <- Postorder(RandomTree(25, root = TRUE))
  for (i in seq_len(12) + 1L) {
    tr[[i]] <- Postorder(TreeSearch::SPR(tr[[i - 1]]))
  }

  write.tree(tr[[3]])
  
  expect_equal(SPRDist(tr[[3]], tr[[11]], method = "DeO"), 8)
  expect_equal(SPRDist(tr[[3]], tr[[11]], method = "rogue")[[1]], 8)
  options(sprH = "ami") # = 8
  # options(sprH = "joint") = 10
  expect_equal(SPRDist(tr[[3]], tr[[11]], method = "confl"), 8)
  
  # TBRDist::USPRDist(tr[[3]], tr[[11]]) = 8
  
  # 
  # 
  # # Simplified example for reproducibility
  # Simplify <- function(tr) {
  #   # Critical to the behaviour: t19, t25, t5
  #   tr |>
  #     DropTip("t17") |> # Reduces by a step
  #     DropTip("t13") |> # Reduces by a step
  #     DropTip("t8") |> # Reduces by a step
  #     DropTip("t2") |> # Reduces by a step
  #     DropTip("t19") |> # Reduces by a step
  #     DropTip(c("t1", "t4", "t6", "t11", "t12", "t14", "t15", "t9", "t24")) # No difference to score
  # }
  # # t3 <- Simplify(tr[[3]])
  # # t11 <- Simplify(tr[[11]])
  # t3 <- Tree("((((t21,(((((t5,t22),t7),t20),t25),(t23,t16))),t18),t10),t3);")
  # t11 <- Tree("((((t21,t18),t10),(((t20,t5),(t7,t22)),(t23,(t16,t25)))),t3);")
  # 
  # TBRDist::USPRDist(t3, t11) # 3
  # expect_equal(SPRDist(t3, t11, method = "DeO"), 3)
  # expect_equal(SPRDist(t3, t11, method = "confl")[[1]], 3)
  # 
  # 
  # 
  # options(debugSPR = T)
  # # Simplified example for reproducibility
  # 
  # Simplify <- function(tr) {
  #   DropTip(tr, paste0("t", c(1:11, 13:14, 16:17, 20, 22, 25)))
  # }
  # # t3 <- Simplify(tr[[3]])
  # # t11 <- Simplify(tr[[11]])
  # write.tree(t3)
  # write.tree(t11)
  # t3 <- Tree( "(((t21,(((t15,t12),t23),t24)),t18),t19);")
  # t11 <- Tree( "((t21,t18),(((t12,(t19,t15)),t23),t24));")
  # for (lab in TipLabels(Simplify(tr[[3]]))) {
  #   d <- SPRDist(DropTip(t3, lab), DropTip(t11, lab), method = "DeO")
  #   c <- SPRDist(DropTip(t3, lab), DropTip(t11, lab), method = "conf")
  #   if (c - d == 3) {
  #     stop(c, ", ", d, ": ", lab)
  #   }
  # }
  # 
  # deO <- SPRDist(t3, t11, method = "DeO")
  # options(sprH = "ami")
  # conf <- SPRDist(t3, t11, method = "confl")
  # TBRDist::USPRDist(t3, t11)
  # message(deO, ", ", conf)
  # expect_equal(conf - deO,
  #              SPRDist(tr[[3]], tr[[11]], method = "confl") -
  #                SPRDist(tr[[3]], tr[[11]], method = "DeO"))
  
  

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
  testDist <- SPRDist(tr, method = "rogue")
  simDist <- dist(seq_along(tr))

  expect_true(all(testDist <= simDist))

  if (interactive()) {
    plot(testDist ~ jitter(simDist), asp = 1, frame.plot = F, col = 3)
    abline(0, 1)
    points(phanDist ~ jitter(simDist), pch = 3, col = 2)
  }
  # bestDist <- as.dist(pmin(as.matrix(testDist), as.matrix(SPRDist(rev(tr)))[rev(seq_along(tr)), rev(seq_along(tr))]))
  bestDist <- testDist # assert symmetry

  overShot <- as.matrix(testDist) > as.matrix(simDist)
  which(overShot, arr.ind = TRUE)
  
  underShot <- as.matrix(testDist) < as.matrix(phanDist)
  which(underShot, arr.ind = TRUE)
  
  if (interactive() && nTip == 25 && nSPR == 12) {
    #  trueDist <- TBRDist::USPRDist(tr)
    trueDist <- readRDS("true-25tip-12spr.Rds")

    par(mfrow = c(1, 2))
    distRange <- c(simDist - phanDist, simDist - bestDist)
    hist(distRange, col = NA, border = NA)
    hist(simDist - phanDist, add = TRUE, col = 2)
    hist(simDist - bestDist, add = TRUE, col = "#88ee4488")
    
    plot(simDist, simDist, type = "n", asp = 1, ylim = range(distRange),
         xlab = "Number of SPR moves", frame.plot = FALSE)
    abline(0, 0, col = 4)
    jd <- jitter(simDist)
    #points(jd, trueDist, pch = 7, col = 3)
    #points(jd, phanDist, pch = 1)
    #points(jd, bestDist, pch = 3, col = 2)
    points(jd, phanDist - trueDist, pch = 5, col = 2)
    points(jd, bestDist - trueDist, pch = 4, col = 3)
    legend("bottomright", c("Phangorn", "Rogue"), bty = "n",
           pch = 5:4, col = 2:3)
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
