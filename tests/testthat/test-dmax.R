dmaxFns <- list(
  CID = list(ClusteringInfoDistance, MutualClusteringInfo, ClusteringEntropy),
  DPI = list(DifferentPhylogeneticInfo, SharedPhylogeneticInfo, SplitwiseInfo),
  MSID = list(MatchingSplitInfoDistance, MatchingSplitInfo, SplitwiseInfo)
)

# Independent per-pair calculation, on explicitly pruned trees
.ExpectedDMax <- function(fns, tree1, tree2) {
  shared <- intersect(TipLabels(tree1), TipLabels(tree2))
  tree1 <- KeepTip(tree1, shared)
  tree2 <- KeepTip(tree2, shared)
  maxInfo <- max(fns[[3]](tree1), fns[[3]](tree2))
  if (maxInfo == 0) 0 else 1 - fns[[2]](tree1, tree2) / maxInfo
}

.ExpectedCross <- function(fns, trees1, trees2) {
  outer(seq_along(trees1), seq_along(trees2), Vectorize(function(i, j) {
    .ExpectedDMax(fns, trees1[[i]], trees2[[j]])
  }))
}

.ExpectedAll <- function(fns, trees) {
  ret <- .ExpectedCross(fns, trees, trees)
  ret[lower.tri(ret)]
}

set.seed(5108)
dmaxTrees <- c(
  bal = BalancedTree(9),
  pec = PectinateTree(9),
  part = CollapseNode(BalancedTree(9), 12:13),
  rnd1 = RandomTree(9, root = TRUE),
  rnd2 = CollapseNode(RandomTree(9, root = TRUE), 14),
  star = StarTree(9)
)
otherTrees <- c(
  o1 = RandomTree(9, root = TRUE),
  o2 = CollapseNode(PectinateTree(9), 15),
  o3 = BalancedTree(9)
)

test_that("normalize = 'dmax' matches per-pair calculation: shared leaves", {
  for (fns in dmaxFns) {
    Dist <- fns[[1]]

    allPairs <- Dist(dmaxTrees, normalize = "dmax")
    expect_s3_class(allPairs, "dist")
    expect_equal(attr(allPairs, "Labels"), names(dmaxTrees))
    expect_equal(as.vector(allPairs), .ExpectedAll(fns, dmaxTrees))

    cross <- Dist(dmaxTrees, otherTrees, normalize = "dmax")
    expect_equal(dimnames(cross), list(names(dmaxTrees), names(otherTrees)))
    expect_equal(unname(cross), .ExpectedCross(fns, dmaxTrees, otherTrees))

    oneMany <- Dist(otherTrees[[1]], dmaxTrees, normalize = "dmax")
    expect_equal(names(oneMany), names(dmaxTrees))
    expect_equal(unname(oneMany),
                 .ExpectedCross(fns, otherTrees[1], dmaxTrees)[1, ])
    manyOne <- Dist(dmaxTrees, otherTrees[[2]], normalize = "dmax")
    expect_equal(names(manyOne), names(dmaxTrees))
    expect_equal(unname(manyOne),
                 .ExpectedCross(fns, dmaxTrees, otherTrees[2])[, 1])

    expect_equal(Dist(dmaxTrees[[1]], dmaxTrees[[2]], normalize = "dmax"),
                 .ExpectedDMax(fns, dmaxTrees[[1]], dmaxTrees[[2]]))

    matched <- Dist(dmaxTrees, normalize = "dmax", reportMatching = TRUE)
    expect_equal(as.vector(matched), as.vector(allPairs))
    pair <- Dist(dmaxTrees[[1]], dmaxTrees[[3]], normalize = "dmax",
                 reportMatching = TRUE)
    expect_false(is.null(attr(pair, "matching")))
    expect_equal(as.vector(pair),
                 .ExpectedDMax(fns, dmaxTrees[[1]], dmaxTrees[[3]]))

    # Identical trees are zero apart; any pair lies in [0, 1]
    expect_equal(Dist(dmaxTrees[[4]], dmaxTrees[[4]], normalize = "dmax"), 0)
    expect_true(all(allPairs >= 0 & allPairs <= 1))
  }
})

test_that("normalize = 'dmax' differs from normalize = TRUE", {
  # The two agree only when both trees carry equal information
  expect_false(isTRUE(all.equal(
    ClusteringInfoDistance(dmaxTrees, normalize = "dmax"),
    ClusteringInfoDistance(dmaxTrees, normalize = TRUE)
  )))
})

test_that("normalize = 'dmax' matches per-pair calculation: differing leaves", {
  trees <- c(
    a = BalancedTree(9),
    b = DropTip(PectinateTree(10), "t10"),
    c = DropTip(RandomTree(11, root = TRUE), c("t1", "t11")),
    d = CollapseNode(DropTip(BalancedTree(10), "t3"), 13)
  )
  others <- c(
    e = DropTip(RandomTree(10, root = TRUE), "t2"),
    f = RandomTree(8, root = TRUE)
  )
  for (fns in dmaxFns) {
    Dist <- fns[[1]]
    expect_equal(as.vector(Dist(trees, normalize = "dmax")),
                 .ExpectedAll(fns, trees))
    expect_equal(unname(Dist(trees, others, normalize = "dmax")),
                 .ExpectedCross(fns, trees, others))
    expect_equal(unname(Dist(trees, others[[1]], normalize = "dmax")),
                 .ExpectedCross(fns, trees, others[1])[, 1])
    expect_equal(unname(Dist(others[[2]], trees, normalize = "dmax")),
                 .ExpectedCross(fns, others[2], trees)[1, ])
    expect_equal(Dist(trees[[1]], others[[2]], normalize = "dmax"),
                 .ExpectedDMax(fns, trees[[1]], others[[2]]))
  }
})

test_that("normalize = 'dmax' is zero between trees with no information", {
  for (fns in dmaxFns) {
    expect_equal(fns[[1]](StarTree(8), StarTree(8), normalize = "dmax"), 0)
    expect_equal(
      as.vector(fns[[1]](c(StarTree(8), StarTree(8)), normalize = "dmax")),
      0
    )
  }
})

test_that("Similarities reject normalize = 'dmax'", {
  for (fns in dmaxFns) {
    Similarity <- fns[[2]]
    expect_error(Similarity(dmaxTrees, normalize = "dmax"), "information-based")
    expect_error(Similarity(dmaxTrees[[1]], dmaxTrees[[2]], normalize = "dmax"),
                 "information-based")
  }
  expect_error(ClusteringInfoDistance(dmaxTrees, normalize = "Dmax"),
               "\"Dmax\".*not available")
})
