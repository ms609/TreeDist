library("TreeTools", quietly = TRUE)

test_that("NNIDist() handles exceptions", {
  expect_error(NNIDist(list(PectinateTree(7), PectinateTree(8))),
               "trees must contain the same number of leaves")
  expect_error(NNIDist(list(PectinateTree(1:8), PectinateTree(8))),
               "trees must bear identical labels")
  expect_error(NNIDist(list(PectinateTree(1:8), 
                            PectinateTree(as.character(1:8)))),
               "trees must bear identical labels")
  expect_error(cpp_nni_distance(
    PectinateTree(40000)$edge, # Will fail before not being postorder is problem
    BalancedTree(40000)$edge, 40000), "so many tips")
  
  expect_error(NNIDist(BalancedTree(5), RootOnNode(BalancedTree(5), 1)))
  
})

test_that("NNIDist() at NNI_MAX_TIPS", {
  maxTips <- 32768
  more <- maxTips + 1
  expect_error(.NNIDistSingle(PectinateTree(more), BalancedTree(more), more),
               "so many tips")
  goingQuickly <- TRUE
  skip_if(goingQuickly)
  
  heapTips <- 16384 + 1
  skip_if_not_installed("testthat", "3.2.2")
  expect_no_error(.NNIDistSingle(PectinateTree(heapTips),
                                 BalancedTree(heapTips), heapTips))
  
  n <- .NNIDistSingle(PectinateTree(maxTips), BalancedTree(maxTips),
                           maxTips)
  expect_gt(n[["best_upper"]], n[["best_lower"]])
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
  
  Fack <- function(n) ((n - 2) * ceiling(log2(n))) + n
                                         
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
  
  allMatched <- c(lower = 0L, best_lower = 0L, tight_upper = 0L,
                  best_upper = 0L, loose_upper = 0L, fack_upper = 0L,
                  li_upper = 0L)
  oneUnmatched <- c(lower = 1L, best_lower = 1L, tight_upper = 1L,
                    best_upper = 1L, loose_upper = 2L, fack_upper = 2L,
                    li_upper = Li(1))
  fiveUnmatched <- c(lower = 5L, best_lower = 10L, tight_upper = 10L,
                     best_upper = 10L, loose_upper = 18L, fack_upper = 18L,
                     li_upper = Li(5))
  
  Test <- function(tree, expect) {
    expectation <- rep(NA_integer_, 7L)
    names(expectation) <- c("lower", "best_lower", "tight_upper", "best_upper",
                            "loose_upper", "fack_upper", "li_upper")
    expectation[names(expect)] <- expect
    if (is.na(expectation["best_lower"]) && !is.na(expect["tight_upper"])) {
      expectation[c("best_lower", "best_upper")] <- expect["tight_upper"]
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
       c(lower = 3L, tight_upper = 5L, fack_upper = 8L, li_upper = Li(3)))
  
  # One region of four unmatched edges
  Test(ape::read.tree(text="(((a, e), (f, d)), ((b, c), (g, h)));"),
       c(lower = 4L, tight_upper = 7L, fack_upper = 14L, li_upper = Li(4)))
  
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
  FackMax <- function(n) n * ceiling(log2(n)) + n - (2 * ceiling(log2(n)))
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
