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

test_that("SPR shortcuts okay - exhaustive", {
  skip_if_not(getOption("slowMode", FALSE))
  skip_if_not_installed("TBRDist")
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
  skip_if_not(getOption("slowMode", FALSE))
  skip_if_not_installed("TBRDist")
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
  skip_if_not(getOption("slowMode", FALSE))
  
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
  
  # Compares against TBRDist's exact uSPR distance, when available (Suggests);
  # the fixed-value expectations below run regardless.
  haveTBRDist <- requireNamespace("TBRDist", quietly = TRUE)
  expect_exact <- function(x, y, method = "rogue") {
    if (!haveTBRDist) {
      return(invisible())
    }
    if (is.character(x)) x <- Tree(x)
    if (is.character(y)) y <- Tree(y)
    expect_equal(
      SPRDist(x, y, method = method)[[1]],
      TBRDist::USPRDist(x, y)
    )
  }
  
  Tree <- function(txt) ape::read.tree(text = txt)

  expect_equal(
    .SPRRogue(
      ape::read.tree(text = "((a, b), (c, d));"),
      ape::read.tree(text = "((a, c), (b, d));")
    )[[1]],
    1L
  )
  expect_equal(
    .SPRRogue(PectinateTree(letters[1:26]), PectinateTree(letters[c(2:26, 1)]))[[1]],
    1L
  )
  
  # Looks simple, but requires ties to be broken suitably
  expect_exact("(a,(d,(b,(c,X))));", "(a,((b,c),(X,d)));", "rogue") # distance = 1
  
  expect_exact("((((b,c),d),e),a);", "(a,(b,((e,c),d)));", "rogue")
  if (haveTBRDist) {
    expect_failure(
      expect_exact("((((b,c),d),e),a);", "(a,(b,((e,c),d)));", "deo")
    )
  }
  
  # Passes with ami, joint, vi
  # Fails with viNorm - should tiebreaker be non-normalized?
  expect_exact("(((t21,(((t15,t12),t23),t24)),t18),t19);",
               "((t21,t18),(((t12,(t19,t15)),t23),t24));", "rogue")
  
  # Example AZ33: IJK and DEF are schlepped
  # Passes with joint, ami, viNorm
  # Fails with vi (5, 6 with tiebreaker), conf (7)
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
  expect_equal(.SPRRogue(lockedMid1, lockedMid2)[[1]], 5)
})

test_that(".BoundedCID() matches ClusteringInfoDist()", {
  library("TreeTools", quietly = TRUE)
  set.seed(174)
  for (nTip in c(3, 4, 9, 25)) {
    tree1 <- RandomTree(nTip, root = TRUE)
    tree2 <- RandomTree(nTip, root = TRUE)
    exact <- ClusteringInfoDist(tree1, tree2)
    # The bounded fill sums its terms in a different order, so agrees with
    # ClusteringInfoDist() only to within rounding error
    expect_equal(.BoundedCID(tree1, tree2), exact, tolerance = 1e-12)
    expect_equal(.BoundedCID(tree1, tree2, exact), exact, tolerance = 1e-12)
    # Within numerical noise of the cutoff: must be exact, to preserve ties
    expect_equal(.BoundedCID(tree1, tree2, exact - 1e-12), exact,
                 tolerance = 1e-12)
    expect_identical(.BoundedCID(tree1, tree2, exact - 1e-12),
                     .BoundedCID(tree1, tree2))
    if (exact > 0) {
      bounded <- .BoundedCID(tree1, tree2, exact / 2)
      expect_gt(bounded, exact / 2)
      expect_lte(bounded, exact + sqrt(.Machine$double.eps))
    }
  }
  
  # Identical trees
  tree <- BalancedTree(12)
  expect_equal(.BoundedCID(tree, tree, 0), 0)
  expect_equal(spr_bounded_cid(Postorder(tree), Postorder(tree), 0), 0)
})

test_that("spr_bounded_cid() matches .BoundedCID()", {
  library("TreeTools", quietly = TRUE)
  set.seed(1174)
  for (nTip in c(3, 4, 9, 25)) {
    tree1 <- Postorder(RandomTree(nTip, root = TRUE))
    tree2 <- Postorder(RenumberTips(RandomTree(nTip, root = TRUE), tree1))
    exact <- .BoundedCID(tree1, tree2)
    for (cutoff in c(Inf, exact, exact - 1e-12, exact / 2, 0)) {
      expect_equal(spr_bounded_cid(tree1, tree2, cutoff),
                   .BoundedCID(tree1, tree2, cutoff))
    }
  }
})

test_that("C++ Rogue scorer matches R implementation", {
  library("TreeTools", quietly = TRUE)
  Rogue <- function(tree1, tree2, ...) {
    opt <- options(...)
    on.exit(options(opt))
    SPRDist(tree1, tree2, method = "Rogue")
  }
  # Dropping one leaf reconciles the trees: abandons scoring early
  tree1 <- BalancedTree(8)
  tree2 <- AddTip(DropTip(tree1, "t1"), 11, "t1")
  expect_equal(spr_rogue_scores(Postorder(tree1),
                                Postorder(RenumberTips(tree2, tree1)), 2L),
               list(scores = c(-Inf, double(7)), pair_scores = NULL))
  for (depth in 1:2) {
    expect_equal(Rogue(tree1, tree2, sprDepth = depth, sprShortcut = 0), 1)
    expect_equal(Rogue(tree1, tree2, sprDepth = depth, sprShortcut = 0,
                       sprCpp = FALSE), 1)
  }
  
  set.seed(2174)
  for (i in 1:12) {
    nTip <- 4 + i
    tree1 <- RandomTree(nTip, root = TRUE)
    tree2 <- RandomTree(nTip, root = TRUE)
    for (depth in 1:2) {
      expect_identical(Rogue(tree1, tree2, sprDepth = depth),
                       Rogue(tree1, tree2, sprDepth = depth, sprCpp = FALSE))
    }
  }
})

test_that("Single-leaf candidates scored from parent splits", {
  library("TreeTools", quietly = TRUE)
  # Scores must agree with those from rebuilding the reduced trees, to within
  # numerical error, with the same best candidates; candidates that cannot
  # win may receive different lower bounds
  Agree <- function(x, y) {
    expect_equal(is.finite(x), is.finite(y))
    fin <- is.finite(x)
    if (!any(fin)) return()
    best <- min(x[fin])
    expect_equal(min(y[fin]), best, tolerance = 1e-12)
    tol <- 1e-9 * max(1, best)
    expect_equal(which(fin & x <= best + tol), which(fin & y <= best + tol))
    exact <- fin & x <= best + 1e-6 & y <= best + 1e-6
    expect_equal(x[exact], y[exact], tolerance = 1e-12)
  }
  Compare <- function(tree1, tree2) {
    reduced <- ReduceTrees(tree1, tree2)
    if (is.null(reduced) || NTip(reduced[[1]]) < 8) return(FALSE)
    for (depth in 1:2) {
      fast <- spr_rogue_scores(reduced[[1]], reduced[[2]], depth)
      slow <- spr_rogue_scores(reduced[[1]], reduced[[2]], depth, FALSE)
      Agree(fast[["scores"]], slow[["scores"]])
      Agree(fast[["pair_scores"]], slow[["pair_scores"]])
    }
    TRUE
  }

  # Removing a leaf leaves a chain of leaves common to both trees, which
  # reduce_trees() would shorten
  chain <- ape::read.tree(text = c(
    "(t1,((((((t2,t11),t10),t15),t13),((((((t3,t8),t7),t14),t4),t6),t5)),(t9,t12)));",
    "(t1,(((t8,((t15,((t2,t11),t10)),t13)),(((((t3,(t7,t9)),t14),t4),t6),t5)),t12));"
  ))
  expect_true(Compare(chain[[1]], chain[[2]]))

  Moved <- function(tree, nMoves) {
    for (i in seq_len(nMoves)) {
      leaf <- sample(tree[["tip.label"]], 1)
      pruned <- DropTip(tree, leaf)
      tree <- AddTip(pruned, sample(pruned[["edge"]][, 2], 1), leaf)
    }
    tree
  }
  set.seed(4174)
  nScored <- 0
  for (i in 1:60) {
    nTip <- sample(8:40, 1)
    tree1 <- RandomTree(nTip, root = TRUE)
    tree2 <- if (i %% 2) RandomTree(nTip, root = TRUE) else
      Moved(tree1, sample(2:5, 1))
    nScored <- nScored + Compare(tree1, tree2)
  }
  expect_gt(nScored, 40)

  # Trees that share a split use the full rebuild for every candidate
  tree1 <- Postorder(RootTree(BalancedTree(12), 1))
  tree2 <- Postorder(RootTree(RenumberTips(PectinateTree(12), tree1), 1))
  expect_equal(spr_rogue_scores(tree1, tree2, 1L),
               spr_rogue_scores(tree1, tree2, 1L, FALSE))
})

test_that("cpp_mutual_clustering_bounded() exits early", {
  library("TreeTools", quietly = TRUE)
  MCI <- function(x, y, minScore) {
    cpp_mutual_clustering_bounded(x, as.Splits(y, x), NTip(x), minScore)
  }
  bal <- as.Splits(BalancedTree(16))
  pec <- as.Splits(PectinateTree(16))
  
  # Unequal entropies: exits before scoring any pair of splits
  exact <- MutualClusteringInfo(bal, pec)
  needed <- (ClusteringEntropy(bal) + ClusteringEntropy(pec)) / 2
  bound <- MCI(bal, pec, needed)
  expect_lt(bound, needed)
  expect_gt(bound, exact)
  expect_equal(MCI(bal, pec, -Inf), exact)
  expect_equal(MCI(bal, pec, exact), exact)
  
  # Equal entropies: exits once a row's best match falls short
  pec2 <- as.Splits(PectinateTree(paste0("t", c(9:16, 1:8))), pec)
  exact <- MutualClusteringInfo(pec, pec2)
  needed <- ClusteringEntropy(pec) - 1e-9
  bound <- MCI(pec, pec2, needed)
  expect_lt(bound, needed)
  expect_gt(bound, exact)
  expect_equal(MCI(pec, pec2, exact), exact)
})

test_that("Rogue SPR results unchanged by early exits", {
  library("TreeTools", quietly = TRUE)
  # Fixed trees: RandomTree() output can vary between TreeTools versions
  trees <- ape::read.tree(text = c(
    "(t1,(((t2,t17),(((((t4,t12),t6),t16),t9),t13)),(((t3,((t7,(((t10,t14),t18),t15)),t8)),t5),t11)));",
    "((((((((t1,t14),(t9,t13)),t15),(t5,t17)),(((t2,t3),t8),(t10,t11))),t4),(t6,(t16,t18))),(t7,t12));",
    "(((((t1,((t2,(((t5,((t8,t12),t11)),(t14,t18)),t6)),(t13,t15))),(t9,t17)),t3),t4),((t7,t10),t16));",
    "(((((((((t1,t8),(t13,t14)),t18),t10),(((((t2,t15),t4),(t6,t12)),t16),t11)),t7),t17),t9),(t3,t5));",
    "((((((t1,((((t3,(t7,t18)),(((((t8,t16),t9),t10),t11),t15)),t6),t12)),t17),t14),t5),t13),(t2,t4));",
    "((((((((((t1,t2),t16),((t3,((t6,t15),t18)),t8)),t9),(((t5,t10),t12),t7)),t14),t11),t4),t13),t17);"
  ))
  Rogue <- function(...) {
    opt <- options(...)
    on.exit(options(opt))
    SPRDist(trees[[1]], trees[-1], method = "Rogue")
  }
  # Values computed before early exits were introduced
  expect_equal(Rogue(sprProxy = "C"), c(10, 10, 11, 10, 11))
  expect_equal(Rogue(sprProxy = "C", sprDepth = 2), c(10, 10, 10, 10, 11))
  expect_equal(Rogue(sprProxy = "P"), c(10, 10, 11, 10, 11))
  expect_equal(Rogue(sprProxy = "R"), c(11, 12, 11, 13, 13))
  expect_equal(Rogue(sprCpp = FALSE), c(10, 10, 11, 10, 11))
  expect_equal(Rogue(sprCpp = FALSE, sprDepth = 2), c(10, 10, 10, 10, 11))
  
  # Abandons a divide-and-conquer candidate after its first half
  pinned1 <- ape::read.tree(text = "(t1,((((t2,t4),t9),(((((((t3,t17),(t5,t7)),t11),t12),t14),t18),(t10,t13))),(t6,((t8,t16),t15))));")
  pinned2 <- ape::read.tree(text = "(((((((((((t1,(t4,t8)),t16),t6),t14),t15),(t2,t12)),t10),((t3,t18),t9)),((t7,t13),t17)),t11),t5);")
  expect_equal(SPRDist(pinned1, pinned2, method = "Rogue"), 9)
  opt <- options(sprCpp = FALSE)
  expect_equal(SPRDist(pinned1, pinned2, method = "Rogue"), 9)
  options(opt)
  skip_if_not_installed("Quartet")
  expect_true(all(Rogue(sprProxy = "Q") %in% 1:15))
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
  expect_equal(SPRDist(PectinateTree(letters[1:26]),
                       PectinateTree(letters[c(2:26, 1)]),
                       method = "deOliv"),
               1L)
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
