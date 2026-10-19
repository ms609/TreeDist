library("TreeTools", quietly = TRUE)

# Reference: prune each pair to its shared leaves, then compare
PrunedOracle <- function(Func, x, y, ...) {
  common <- intersect(TipLabels(x), TipLabels(y))
  Func(KeepTip(x, common), KeepTip(y, common), ...)
}

OracleAll <- function(Func, trees, ...) {
  n <- length(trees)
  pairs <- combn(n, 2)
  apply(pairs, 2, function(ij) {
    PrunedOracle(Func, trees[[ij[[1]]]], trees[[ij[[2]]]], ...)
  })
}

# Pairwise values, in combn() order, from a dist object or a full matrix
AsPairs <- function(x) {
  if (inherits(x, "dist")) as.numeric(x) else x[lower.tri(x)]
}

OracleCross <- function(Func, trees1, trees2, ...) {
  matrix(
    vapply(seq_along(trees2), function(j) {
      vapply(trees1, PrunedOracle, double(1), Func = Func, y = trees2[[j]],
             ...)
    }, double(length(trees1))),
    length(trees1), length(trees2),
    dimnames = list(names(trees1), names(trees2))
  )
}

{ # Reference objects
  set.seed(8231)
  labs <- paste0("t", 1:24)
  RandTree <- function(n) {
    tr <- RandomTree(sample(labs, n), root = TRUE)
    if (n > 8) {
      # Polytomies
      tr <- CollapseNode(tr, n + 1 + sample.int(tr[["Nnode"]] - 1, 2))
    }
    tr
  }
  trees <- structure(lapply(c(24, 20, 16, 12, 9, 6), RandTree),
                     names = letters[1:6], class = "multiPhylo")
  others <- structure(lapply(c(22, 14, 8), RandTree),
                      names = LETTERS[1:3], class = "multiPhylo")
}

test_that("Distances between trees with different leaves match pruned trees", {
  withr::local_options(mc.cores = 2L)
  Jaccard2 <- function(...) {
    JaccardRobinsonFoulds(..., k = 2, allowConflict = FALSE)
  }
  RFSimilarity <- function(...) RobinsonFoulds(..., similarity = TRUE)
  fns <- list(ClusteringInfoDistance, MutualClusteringInfo,
              DifferentPhylogeneticInfo, SharedPhylogeneticInfo,
              MatchingSplitInfoDistance, MatchingSplitInfo,
              MatchingSplitDistance, InfoRobinsonFoulds, RobinsonFoulds,
              RFSimilarity, NyeSimilarity, Jaccard2)
  for (Func in fns) {
    # MatchingSplitDistance has no default normalization
    norms <- if (identical(Func, MatchingSplitDistance)) FALSE else c(FALSE, TRUE)
    for (normalize in norms) {
      expect_equal(AsPairs(Func(trees, normalize = normalize)),
                   OracleAll(Func, trees, normalize = normalize))
      expect_equal(Func(trees, others, normalize = normalize),
                   OracleCross(Func, trees, others, normalize = normalize))
      expect_equal(Func(trees, others[[1]], normalize = normalize),
                   OracleCross(Func, trees, others[1],
                               normalize = normalize)[, 1])
      expect_equal(Func(others[[2]], trees, normalize = normalize),
                   OracleCross(Func, trees, others[2],
                               normalize = normalize)[, 1])
    }
  }

  # Normalizing against a function of each tree's information
  expect_equal(MutualClusteringInfo(trees, others, normalize = pmin),
               OracleCross(MutualClusteringInfo, trees, others,
                           normalize = pmin))
  expect_equal(AsPairs(SharedPhylogeneticInfo(trees, normalize = pmax)),
               OracleAll(SharedPhylogeneticInfo, trees, normalize = pmax))
})

test_that("RobinsonFoulds() detects trees whose leaves are a subset", {
  # The first tree in each list holds every leaf
  full <- trees[c(1, 1)]
  part <- trees[c(1, 3)]
  expect_equal(RobinsonFoulds(full, part),
               OracleCross(RobinsonFoulds, full, part))
  expect_equal(RobinsonFoulds(part, full, similarity = TRUE),
               OracleCross(RobinsonFoulds, part, full, similarity = TRUE))
  expect_equal(AsPairs(RobinsonFoulds(part)),
               OracleAll(RobinsonFoulds, part))
  expect_null(.SameLeafClusterTables(part))
  expect_length(.SameLeafClusterTables(full), 2)
  expect_length(.SameLeafClusterTables(trees[[1]]), 1)
  expect_null(.SameLeafClusterTables(trees[c(3, 1)]))
})

test_that("Pruned pairs: results are independent of thread count", {
  withr::local_options(mc.cores = 1L)
  serial <- TreeDistance(trees)
  withr::local_options(mc.cores = 2L)
  expect_equal(TreeDistance(trees), serial)
})

test_that("Pruned pairs: Splits input", {
  splits <- lapply(trees, as.Splits)
  expect_equal(
    AsPairs(CalculateTreeDistance(MutualClusteringInfoSplits, splits)),
    OracleAll(MutualClusteringInfo, trees)
  )
  expect_equal(
    unname(CalculateTreeDistance(MutualClusteringInfoSplits, splits[1:2],
                                 splits[3:5])),
    unname(OracleCross(MutualClusteringInfo, trees[1:2], trees[3:5]))
  )
})

test_that("Pruned pairs: splits that coincide after pruning are merged", {
  # Dropping X makes ((A, B), X) and (A, B) identical
  tr1 <- ape::read.tree(text = "((((A, B), X), C), (D, (E, F)));")
  tr2 <- ape::read.tree(text = "(((A, B), C), (D, (E, (F, Y))));")
  tr3 <- ape::read.tree(text = "(((A, C), B), ((D, E), F));")
  trees <- c(tr1, tr2, tr3)
  expect_equal(AsPairs(MatchingSplitDistance(trees)),
               OracleAll(MatchingSplitDistance, trees))
  for (Func in list(ClusteringInfoDistance, SharedPhylogeneticInfo,
                    RobinsonFoulds)) {
    expect_equal(AsPairs(Func(trees)), OracleAll(Func, trees))
    expect_equal(AsPairs(Func(trees, normalize = TRUE)),
                 OracleAll(Func, trees, normalize = TRUE))
  }
  expect_equal(.PrunedPairs(tr1, tr2, Value = NSplits)[-1],
               list(info1 = 3, info2 = 3))
  expect_equal(RobinsonFoulds(tr1, tr2), 0)
})

test_that("Pruned pairs: fewer than four shared leaves", {
  tr1 <- BalancedTree(c("a", "b", "c", "x1", "x2", "x3"))
  tr2 <- PectinateTree(c("a", "b", "c", "y1", "y2", "y3"))
  tr3 <- BalancedTree(c("x1", "x2", "y1", "y2", "z"))
  trees <- c(tr1, tr2, tr3)
  expect_equal(AsPairs(MutualClusteringInfo(trees)), c(0, 0, 0))
  expect_equal(AsPairs(ClusteringInfoDistance(trees)), c(0, 0, 0))
  # No information to normalize against
  expect_equal(AsPairs(TreeDistance(trees)), rep(NaN, 3))
  expect_equal(
    .PrunedPairs(trees, Func = SharedPhylogeneticInfoSplits,
                 Value = SplitwiseInfo),
    list(score = c(0, 0, 0), info1 = c(0, 0, 0), info2 = c(0, 0, 0))
  )
  expect_equal(.PrunedPairs(trees, Value = ClusteringEntropy)[["info1"]],
               c(0, 0, 0))
})

test_that("Pruned pairs: trees with more than 64 leaves", {
  skip_on_cran()
  set.seed(9214)
  bigLabs <- paste0("t", 1:140)
  big <- structure(lapply(c(140, 100, 75), function(n) {
    RandomTree(sample(bigLabs, n), root = TRUE)
  }), class = "multiPhylo")
  for (Func in list(ClusteringInfoDistance, PhylogeneticInfoDistance,
                    RobinsonFoulds)) {
    expect_equal(AsPairs(Func(big, normalize = TRUE)),
                 OracleAll(Func, big, normalize = TRUE))
  }
})

test_that("Pruned pairs: empty input", {
  expect_equal(.PrunedPairs(trees[1], Func = MutualClusteringInfoSplits),
               list(score = double(0), info1 = double(0), info2 = double(0)))
})

test_that("R fallbacks for functions without a C++ implementation", {
  Leaves <- function(tree) NTip(tree)
  pairs <- combn(length(trees), 2)
  shared <- apply(pairs, 2, function(ij) {
    length(intersect(TipLabels(trees[[ij[1]]]), TipLabels(trees[[ij[2]]])))
  })
  expect_equal(.MaxValue(trees, NULL, Leaves), shared + shared)

  sharedCross <- outer(seq_along(trees), seq_along(others),
                       Vectorize(function(i, j) {
                         length(intersect(TipLabels(trees[[i]]),
                                          TipLabels(others[[j]])))
                       }))
  expect_equal(.MaxValue(trees, others, Leaves),
               as.numeric(sharedCross + sharedCross))

  expect_equal(
    NormalizeInfo(matrix(1, 6, 3), trees, others, InfoInTree = Leaves),
    1 / (sharedCross + sharedCross)
  )
  # Additional arguments to InfoInTree
  expect_equal(
    NormalizeInfo(matrix(1, 6, 3), trees, others, InfoInTree = SplitwiseInfo,
                  sum = TRUE),
    unname(1 / OracleCross(function(x, y) SplitwiseInfo(x) + SplitwiseInfo(y),
                           trees, others))
  )
  # Numeric `how` needs no pruning
  expect_equal(NormalizeInfo(matrix(1, 6, 3), trees, others, how = 4),
               matrix(0.25, 6, 3))

  # Scoring function without a C++ implementation: rows match tree1
  NSharedSplits <- function(splits1, splits2, nTip, ...) nTip
  expect_equal(CalculateTreeDistance(NSharedSplits, trees, others),
               sharedCross, ignore_attr = TRUE)
  expect_equal(dimnames(CalculateTreeDistance(NSharedSplits, trees, others)),
               list(names(trees), names(others)))
})
