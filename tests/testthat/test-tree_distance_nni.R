library("TreeTools", quietly = TRUE)

test_that("NNIDist() handles exceptions", {
  expect_error(NNIDist(list(PectinateTree(7), PectinateTree(8))),
               "trees must contain the same number of leaves")
  expect_error(NNIDist(list(PectinateTree(1:8), PectinateTree(8))),
               "trees must bear identical labels")
  expect_error(NNIDist(list(PectinateTree(1:8),
                            PectinateTree(as.character(1:8)))),
               "trees must bear identical labels")
  expect_error(NNIDist(BalancedTree(5), RootOnNode(BalancedTree(5), 1)))
})

test_that("NNIDist() rejects too many tips (TreeTools >= 2.3.0)", {
  skip_if(TreeDist:::cpp_sl_max_tips() <= 2048L,
          "Compiled against TreeTools < 2.3.0 (SL_MAX_TIPS = 2048)")
  expect_error(NNIDist(PectinateTree(40000), BalancedTree(40000)),
               "not yet supported")
  expect_error(.NNIDistSingle(PectinateTree(32769), BalancedTree(32769), 32769L),
               "not yet supported")
})

test_that("NNIDist() rejects too many tips (TreeTools < 2.3.0)", {
  skip_if(TreeDist:::cpp_sl_max_tips() > 2048L,
          "Compiled against TreeTools >= 2.3.0 (SL_MAX_TIPS > 2048)")
  expect_error(NNIDist(PectinateTree(40000), BalancedTree(40000)),
               "exceed the compiled limit")
  expect_error(.NNIDistSingle(PectinateTree(32769), BalancedTree(32769), 32769L),
               "exceed the compiled limit")
})

test_that("NNIDist() at NNI_MAX_TIPS", {
  maxTips <- 32768
  goingQuickly <- TRUE
  skip_if(goingQuickly)
  
  heapTips <- 16384 + 1
  skip_if_not_installed("testthat", "3.2.2")
  expect_no_error(.NNIDistSingle(PectinateTree(heapTips),
                                 BalancedTree(heapTips), heapTips))
  
  n <- .NNIDistSingle(PectinateTree(maxTips), BalancedTree(maxTips),
                           maxTips)
  expect_gt(n[["best_upper"]], n[["lower"]])
  if (!is.na(n[["tight_upper"]])) {
    expect_gte(n[["tight_upper"]], n[["best_upper"]])
  }
  expect_gte(n[["loose_upper"]], n[["best_upper"]])
  expect_gte(n[["fack_upper"]], n[["best_upper"]])
  expect_gte(n[["li_upper"]], n[["best_upper"]])
})

test_that("Simple NNI approximations", {
  nTip <- 6L
  tree1 <- BalancedTree(nTip)
  tree2 <- PectinateTree(nTip)
  edge1 <- Postorder(tree1$edge)
  edge2 <- Postorder(tree2$edge)
  
  # Fack et al. 2002, theorem 3; N = internal nodes = n_edges + 1
  Fack <- function(n_edges) {
    N <- n_edges + 1
    lgN <- ceiling(log2(N))
    N * lgN - 2 ^ lgN + 1 + 2 * (N - ceiling(log2(N + 1)))
  }
                                         
  Sorting <- function(n) {
    lc <- ceiling(log2(n))
    n * lc - 2 ^ lc + 1
  }
  DegDist <- function(n) {
    nif <- ceiling(log2(n / 3))
    tif <- 2 ^ nif
    tl = n - tif
    mbn <- nif + ceiling(log2(tl / 2)) + 1
    n - 2 - mbn
  }
  Li <- function(n_edges) Sorting(n_edges + 3) + (2 * DegDist(n_edges + 3))
  
  allMatched <- c(lower = 0L, best_upper = 0L, tight_upper = 0L,
                  loose_upper = 0L, fack_upper = 0L, li_upper = 0L)
  oneUnmatched <- c(lower = 1L, best_upper = 1L, tight_upper = 1L,
                    loose_upper = 1L, fack_upper = Fack(1), li_upper = Li(1))
  fiveUnmatched <- c(lower = 5L, best_upper = 10L, tight_upper = 10L,
                     loose_upper = 17L, fack_upper = Fack(5), li_upper = Li(5))
  
  Test <- function(tree, expect) {
    expectation <- rep(NA_integer_, 6L)
    names(expectation) <- c("lower", "best_upper", "tight_upper",
                            "loose_upper", "fack_upper", "li_upper")
    expectation[names(expect)] <- expect
    if (is.na(expectation["best_upper"]) && !is.na(expect["tight_upper"])) {
      expectation["best_upper"] <- expect["tight_upper"]
    }
    if (is.na(expect["loose_upper"])) {
      expectation["loose_upper"] <- min(expect[c("fack_upper", "li_upper")])
    }
    
    expect_equal(NNIDist(tree1, tree), expectation)
    for (i in c(2L, 3L, 4L, 6L)) {
      tree1i <- RootOnNode(tree1, i)
      j <- 0
      for (t2 in unique(lapply(1:9, RootOnNode, tree = tree))) {
        expect_equal(NNIDist(tree1i, t2), expectation)
      }
    }
  }
  
  expect_equal(NNIDist(BalancedTree(2), PectinateTree(2)), allMatched)
  
  expect_equal(cpp_nni_distance(edge1, edge2, NTip(tree1)), oneUnmatched)
  Test(PectinateTree(nTip), oneUnmatched)

  # Identical trees
  tree1 <- Postorder(read.tree(text = "(((a, b), (c, d)), ((e, f), (g, h)));"))
  tree2 <- Postorder(read.tree(text = "(((a, b), (d, c)), ((h, g), (f, e)));"))
  Test(tree1, allMatched)
  Test(tree2, allMatched)
  
  # Tree names
  output <- NNIDist(list(bal = tree1, pec = tree2), 
                    as.phylo(0:2, tipLabels = letters[1:8]))
  expect_equal(rownames(output), c("bal", "pec"))
  
  # Only root edge is different
  Test(Postorder(ape::read.tree(text="(((a, b), (e, f)), ((c, d), (g, h)));")),
       oneUnmatched)
  
  # Two separate regions of difference one
  Test(read.tree(text="((((a, b), c), d), (e, (f, (g, h))));"),
       oneUnmatched * 2)
  
  # One region of three unmatched edges
  Test(read.tree(text="(((a, e), (c, d)), ((b, f), (g, h)));"),
       c(lower = 3L, tight_upper = 5L, fack_upper = Fack(3), li_upper = Li(3)))
  
  # One region of four unmatched edges
  Test(ape::read.tree(text="(((a, e), (f, d)), ((b, c), (g, h)));"),
       c(lower = 4L, tight_upper = 7L, fack_upper = Fack(4), li_upper = Li(4)))
  
  # One region of five unmatched edges
  Test(ape::read.tree(text="(((a, e), (f, d)), ((b, g), (c, h)));"),
       fiveUnmatched)
  
  # Trees with different leaves at root.
  tree1 <- PectinateTree(1:8) # used in Test()
  Test(ape::read.tree(text = "(3, ((5, 6), (7, (1, (2, (4, 8))))));"),
       fiveUnmatched)
  
  # Too different for tight upper bound
  expect_true(is.na(NNIDist(BalancedTree(100), 
                            PectinateTree(100))["tight_upper"]))
  
  # Large, different trees: check that 64 leaf disagreements don't cause crash
  expect_gt(NNIDist(RandomTree(600), RandomTree(600))["li_upper"], 1)
  
})

test_that("NNIDist() lower bounds are valid", {
  # Exact NNI distance is 2 (ms609/TreeDist issue 220)
  a <- read.tree(text = "(t1,(t2,(t4,t5)),(t3,t6));")
  b <- read.tree(text = "(t6,(t1,(t4,t5)),(t3,t2));")
  expect_equal(NNIDist(a, b)[c("lower")], c(lower = 2L))
})

.ExpectBoundsBracketExact <- function(nTip) {
  # Exact distances by breadth-first search over all trees, which are an NNI
  # apart iff they differ in a single split
  trees <- as.phylo(seq_len(NUnrooted(nTip)) - 1, nTip)
  adjacent <- as.matrix(RobinsonFoulds(trees)) == 2
  exact <- matrix(NA_integer_, length(trees), length(trees))
  for (origin in seq_along(trees)) {
    dist <- rep(NA_integer_, length(trees))
    dist[origin] <- 0L
    frontier <- origin
    step <- 0L
    while (length(frontier)) {
      step <- step + 1L
      reached <- which(colSums(adjacent[frontier, , drop = FALSE]) > 0 &
                         is.na(dist))
      dist[reached] <- step
      frontier <- reached
    }
    exact[origin, ] <- dist
  }
  testthat::expect_equal(max(exact), NNIDiameter(nTip)[, "exact"],
                         ignore_attr = TRUE)

  testthat::expect_equal(
    vapply(trees, .NNIExact, integer(1), tree1 = trees[[nTip]]),
    exact[nTip, ]
  )

  bounds <- lapply(NNIDist(trees, regionTips = nTip), as.matrix)
  exact <- exact[lower.tri(exact)]
  .Bound <- function(field) bounds[[field]][lower.tri(bounds[[field]])]
  testthat::expect_true(all(.Bound("lower") <= exact))
  for (upper in c("region_upper", "tight_upper", "best_upper", "loose_upper",
                  "fack_upper", "li_upper")) {
    testthat::expect_true(all(.Bound(upper) >= exact))
  }
  testthat::expect_true(all(.Bound("region_upper") <= .Bound("best_upper")))
  # Not guaranteed in general, but holds for every pair of <= 7 leaves
  testthat::expect_equal(.Bound("region_upper"), exact)
}

test_that("NNIDist() bounds bracket exact distances: 6 leaves", {
  .ExpectBoundsBracketExact(6)
})

test_that("NNIDist() bounds bracket exact distances: 7 leaves", {
  skip_on_cran() # c. 70 s
  .ExpectBoundsBracketExact(7)
})

test_that("NNIDist() region_upper", {
  bal <- BalancedTree(9)
  pec <- PectinateTree(9)
  expect_false("region_upper" %in% names(NNIDist(bal, pec)))
  expect_equal(NNIDist(bal, pec, regionTips = 9)[-7], NNIDist(bal, pec))
  expect_equal(NNIDist(bal, bal, regionTips = 9)[["region_upper"]], 0L)
  expect_equal(
    NNIDist(PectinateTree(3), PectinateTree(3), regionTips = 9),
    c(lower = 0L, best_upper = 0L, tight_upper = 0L, loose_upper = 0L,
      fack_upper = 0L, li_upper = 0L, region_upper = 0L)
  )

  # Two regions, of 7 and 9 leaves, separated by a shared split
  aSide1 <- "((a1, a2), (a3, a4)), (a5, a6)"
  aSide2 <- "(a1, (a3, (a5, (a2, a4)))), a6"
  bSide1 <- "((b1, b2), (b3, b4)), ((b5, b6), (b7, b8))"
  bSide2 <- "b1, (b3, (b5, (b7, (b2, (b4, (b6, b8))))))"
  .Tree <- function(a, b) read.tree(text = paste0("(", a, ", (", b, "));"))
  tree1 <- .Tree(aSide1, bSide1)
  tree2 <- .Tree(aSide2, bSide2)
  bounds <- NNIDist(tree1, tree2, regionTips = 16)
  expect_equal(bounds[["lower"]], 10L)
  expect_equal(bounds[["best_upper"]], 19L)
  expect_equal(bounds[["region_upper"]], 14L) # == .NNIExact(); c. 8 s
  expect_equal(
    NNIDist(RootTree(tree1, "a1"), RootTree(tree2, "b3"), regionTips = 16),
    bounds
  )
  bSide <- paste0("b", 1:8)
  expect_equal(
    NNIDist(RootTree(tree1, bSide), RootTree(tree2, bSide), regionTips = 16),
    bounds
  )
  expect_equal(NNIDist(tree2, tree1, regionTips = 16), bounds)
  expect_equal(
    NNIDist(tree1, list(tree1, tree2), regionTips = 16)["region_upper", ],
    c(0L, 14L)
  )

  # Regions above the cap contribute best_upper
  expect_equal(
    NNIDist(tree1, tree2, regionTips = 6)[["region_upper"]],
    bounds[["best_upper"]]
  )
  expect_equal(NNIDist(tree1, tree2, regionTips = 7)[["region_upper"]], 17L)
  set.seed(0)
  repeat { # Single region of 20 leaves
    big1 <- RandomTree(20, root = FALSE)
    big2 <- RenumberTips(RandomTree(20, root = FALSE), big1)
    big <- NNIDist(big1, big2, regionTips = 16)
    if (big[["lower"]] == 17) break
  }
  expect_equal(big[["region_upper"]], big[["best_upper"]])

  # Polytomies fall back to best_upper, whichever tree holds them
  poly1 <- .Tree(sub("(a1, a2)", "a1, a2", aSide1, fixed = TRUE), bSide1)
  poly2 <- .Tree(aSide2, sub("(b6, b8)", "b6, b8", bSide2, fixed = TRUE))
  polyBounds <- NNIDist(poly1, poly2, regionTips = 16)
  expect_equal(polyBounds[["region_upper"]], polyBounds[["best_upper"]])

  expect_equal(.NNIExact(PectinateTree(3), PectinateTree(3)), 0L)
  expect_equal(.NNIExact(bal, bal), 0L)
  expect_equal(.NNIExact(bal, pec), 2L)
  expect_true(is.na(.NNIExact(BalancedTree(12), PectinateTree(12), 10)))
  expect_error(.NNIExact(PectinateTree(17), PectinateTree(17)), "16 leaves")
  expect_error(.NNIExact(poly1, poly2), "binary")
})

test_that("NNI with lists of trees", {
  tree1 <- BalancedTree(1:8)
  list1 <- list(tree1, PectinateTree(as.character(1:8)),
                PectinateTree(as.character(c(4:1, 5:8))),
                BalancedTree(c(1:3, 8:4)))
  
  multResult <- NNIDist(tree1, list1)
  expect_equal(NNIDist(tree1, list1[[1]]), multResult[, 1])
  expect_equal(NNIDist(tree1, list1[[2]]), multResult[, 2])
  expect_equal(NNIDist(tree1, list1[[3]]), multResult[, 3])
  expect_equal(NNIDist(tree1, list1[[4]]), multResult[, 4])
  
  expect_equal(NNIDist(tree1, list1), NNIDist(list1, tree1))
  
  # CompareAll
  expect_equal(CompareAll(list1, NNIDist), NNIDist(list1))
  
  expect_equal(
    vapply(NNIDist(list1), function(x) unname(as.matrix(x)[1:4, 4:1]),
           matrix(0,4,4)),
    NNIDist(list1, rev(list1)),
    ignore_attr = TRUE
  )
})

test_that("NNIDiameter() is sane", {
  
  exacts <- NNIDiameter(3:12)
  expect_equal(do.call(rbind, NNIDiameter(lapply(3:12, as.integer))), exacts)
  expect_true(all(exacts[, "min"] <= exacts[, "exact"]))
  expect_true(all(exacts[, "max"] >= exacts[, "exact"]))
  expect_true(is.na(NNIDiameter(13)[, "exact"]))
  expect_true(is.na(NNIDiameter(1)[, "exact"]))
  expect_equal(NNIDiameter(BalancedTree(8))[, "exact"], c(exact = 10L))
  
  FackMin <- function(n) ceiling(0.25 * lfactorial(n) / log(2))
  exacts <- c(0, 0, 0, 1, 3, 5, 7, 10, 12, 15, 18, 21)
  liMaxes <- c(0, 1, 3, 5, 8, 13, 16, 21, 25, 31, 37, 43, 47, 53, 59, 65)
  FackMax <- function(n) {
    lgN <- ceiling(log2(n))
    n * lgN - 2 ^ lgN + 1 + 2 * (n - ceiling(log2(n + 1)))
  }
  n <- 4:8
  expect_equal(NNIDiameter(n), cbind(
    liMin = n - 3L,
    fackMin = FackMin(n - 2L),
    min = pmax(n - 3L, FackMin(4:8 - 2L)),
    exact = exacts[n],
    liMax = liMaxes[n],
    fackMax = FackMax(n - 2L),
    max = pmin(liMaxes[n], FackMax(n - 2L))
  ))

  expect_equal(NNIDiameter(as.phylo(0:1, 6)), NNIDiameter(c(6, 6)))
})
