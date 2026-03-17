## Coverage tests for batch C++ functions, fast paths, and cross-pairs.
##
## These tests target code paths introduced during the OpenMP / sort+merge /
## CostMatrix-pooling optimisation cycle that are not exercised by the
## existing test suite.

library("TreeTools", quietly = TRUE)

# Shared fixtures ----
tips20 <- paste0("t", seq_len(20))
tA <- ape::as.phylo(0:4, tipLabels = tips20)
tB <- ape::as.phylo(5:9, tipLabels = tips20)
trees20 <- ape::as.phylo(0:9, tipLabels = tips20)

tips8 <- paste0("t", seq_len(8))
tA8 <- ape::as.phylo(0:2, tipLabels = tips8)
tB8 <- ape::as.phylo(3:5, tipLabels = tips8)


# DifferentPhylogeneticInfo ----

test_that("DPI all-pairs fast path agrees with per-pair", {
  batch <- DifferentPhylogeneticInfo(trees20)
  m <- as.matrix(batch)
  expect_equal(m[2, 1],
               DifferentPhylogeneticInfo(trees20[[1]], trees20[[2]]),
               tolerance = 1e-10)
  expect_equal(m[5, 3],
               DifferentPhylogeneticInfo(trees20[[3]], trees20[[5]]),
               tolerance = 1e-10)
})

test_that("DPI cross-pairs fast path agrees with per-pair", {
  cross <- DifferentPhylogeneticInfo(tA, tB)
  expect_equal(dim(cross), c(5L, 5L))
  expect_equal(cross[1, 1],
               DifferentPhylogeneticInfo(tA[[1]], tB[[1]]),
               tolerance = 1e-10)
  expect_equal(cross[3, 2],
               DifferentPhylogeneticInfo(tA[[3]], tB[[2]]),
               tolerance = 1e-10)
})


# MatchingSplitInfoDistance ----

test_that("MSID all-pairs fast path agrees with per-pair", {
  batch <- MatchingSplitInfoDistance(trees20)
  m <- as.matrix(batch)
  expect_equal(m[2, 1],
               MatchingSplitInfoDistance(trees20[[1]], trees20[[2]]),
               tolerance = 1e-10)
  expect_equal(m[5, 3],
               MatchingSplitInfoDistance(trees20[[3]], trees20[[5]]),
               tolerance = 1e-10)
})

test_that("MSID cross-pairs fast path agrees with per-pair", {
  cross <- MatchingSplitInfoDistance(tA, tB)
  expect_equal(dim(cross), c(5L, 5L))
  expect_equal(cross[1, 1],
               MatchingSplitInfoDistance(tA[[1]], tB[[1]]),
               tolerance = 1e-10)
  expect_equal(cross[3, 2],
               MatchingSplitInfoDistance(tA[[3]], tB[[2]]),
               tolerance = 1e-10)
})


# InfoRobinsonFoulds (distance mode) ----

test_that("IRF all-pairs fast path agrees with per-pair", {
  batch <- InfoRobinsonFoulds(trees20)
  m <- as.matrix(batch)
  expect_equal(m[2, 1],
               InfoRobinsonFoulds(trees20[[1]], trees20[[2]]),
               tolerance = 1e-10)
  expect_equal(m[5, 3],
               InfoRobinsonFoulds(trees20[[3]], trees20[[5]]),
               tolerance = 1e-10)
})

test_that("IRF cross-pairs fast path agrees with per-pair", {
  cross <- InfoRobinsonFoulds(tA, tB)
  expect_equal(dim(cross), c(5L, 5L))
  expect_equal(cross[1, 1],
               InfoRobinsonFoulds(tA[[1]], tB[[1]]),
               tolerance = 1e-10)
  expect_equal(cross[3, 2],
               InfoRobinsonFoulds(tA[[3]], tB[[2]]),
               tolerance = 1e-10)
})


# ClusteringInfoDistance (all-pairs already tested; add cross-pairs variety) ----

test_that("CID cross-pairs with different-sized collections", {
  tA3 <- ape::as.phylo(0:2, tipLabels = tips20)
  tB7 <- ape::as.phylo(3:9, tipLabels = tips20)
  cross <- ClusteringInfoDistance(tA3, tB7)
  expect_equal(dim(cross), c(3L, 7L))
  expect_equal(cross[2, 5],
               ClusteringInfoDistance(tA3[[2]], tB7[[5]]),
               tolerance = 1e-10)
})


# MSD cross-pairs via CalculateTreeDistance → .SplitDistanceManyMany ----

test_that("MSD cross-pairs batch agrees with per-pair", {
  cross <- MatchingSplitDistance(tA8, tB8)
  expect_equal(dim(cross), c(3L, 3L))
  expect_equal(cross[1, 1],
               MatchingSplitDistance(tA8[[1]], tB8[[1]]),
               tolerance = 1e-10)
  expect_equal(cross[2, 3],
               MatchingSplitDistance(tA8[[2]], tB8[[3]]),
               tolerance = 1e-10)
})


# Nye / Jaccard cross-pairs via CalculateTreeDistance → .SplitDistanceManyMany ----

test_that("NyeSimilarity cross-pairs batch agrees with per-pair", {
  cross <- NyeSimilarity(tA8, tB8)
  expect_equal(dim(cross), c(3L, 3L))
  expect_equal(cross[1, 1],
               NyeSimilarity(tA8[[1]], tB8[[1]]),
               tolerance = 1e-10)
})

test_that("JRF cross-pairs batch with k and allowConflict", {
  # Default k=1, allowConflict=TRUE
  cross1 <- JaccardRobinsonFoulds(tA8, tB8)
  expect_equal(dim(cross1), c(3L, 3L))
  expect_equal(cross1[2, 1],
               JaccardRobinsonFoulds(tA8[[2]], tB8[[1]]),
               tolerance = 1e-10)

  # k = 2
  cross_k2 <- JaccardRobinsonFoulds(tA8, tB8, k = 2)
  expect_equal(cross_k2[1, 1],
               JaccardRobinsonFoulds(tA8[[1]], tB8[[1]], k = 2),
               tolerance = 1e-10)

  # allowConflict = FALSE
  cross_nc <- JaccardRobinsonFoulds(tA8, tB8, allowConflict = FALSE)
  expect_equal(cross_nc[1, 1],
               JaccardRobinsonFoulds(tA8[[1]], tB8[[1]],
                                     allowConflict = FALSE),
               tolerance = 1e-10)

  # k = Inf
  cross_inf <- JaccardRobinsonFoulds(tA8, tB8, k = Inf)
  expect_equal(cross_inf[2, 2],
               JaccardRobinsonFoulds(tA8[[2]], tB8[[2]], k = Inf),
               tolerance = 1e-10)
})


# MSI cross-pairs via CalculateTreeDistance → .SplitDistanceManyMany ----

test_that("MSI cross-pairs batch agrees with per-pair", {
  cross <- MatchingSplitInfo(tA8, tB8)
  expect_equal(dim(cross), c(3L, 3L))
  expect_equal(cross[1, 1],
               MatchingSplitInfo(tA8[[1]], tB8[[1]]),
               tolerance = 1e-10)
})

# SPI cross-pairs via CalculateTreeDistance → .SplitDistanceManyMany ----

test_that("SPI cross-pairs batch agrees with per-pair", {
  cross <- SharedPhylogeneticInfo(tA8, tB8)
  expect_equal(dim(cross), c(3L, 3L))
  expect_equal(cross[1, 1],
               SharedPhylogeneticInfo(tA8[[1]], tB8[[1]]),
               tolerance = 1e-10)
})


# C++ batch entropy/info functions ----

test_that("cpp_clustering_entropy_batch matches R ClusteringEntropy", {
  trees <- ape::as.phylo(0:4, tipLabels = tips20)
  splits_list <- as.Splits(trees)
  nTip <- length(tips20)
  batch <- TreeDist:::cpp_clustering_entropy_batch(splits_list, as.integer(nTip))
  r_ref <- vapply(trees, ClusteringEntropy, double(1))
  expect_equal(batch, unname(r_ref), tolerance = 1e-12)
})

test_that("cpp_splitwise_info_batch matches R SplitwiseInfo", {
  trees <- ape::as.phylo(0:4, tipLabels = tips20)
  splits_list <- as.Splits(trees)
  nTip <- length(tips20)
  batch <- TreeDist:::cpp_splitwise_info_batch(splits_list, as.integer(nTip))
  r_ref <- vapply(trees, SplitwiseInfo, double(1))
  expect_equal(batch, unname(r_ref), tolerance = 1e-8)
})

test_that("cpp_clustering_entropy_batch: polytomy trees", {
  # Trees with a single non-trivial split (polytomy)
  poly1 <- ape::read.tree(text = "((t1,t2,t3),(t4,t5,t6));")
  poly2 <- ape::read.tree(text = "((t1,t2),(t3,t4,t5,t6));")
  trees <- structure(list(poly1, poly2), class = "multiPhylo")
  splits_list <- as.Splits(trees)
  nTip <- 6L
  batch <- TreeDist:::cpp_clustering_entropy_batch(splits_list, nTip)
  r_ref <- vapply(trees, ClusteringEntropy, double(1))
  expect_equal(batch, unname(r_ref), tolerance = 1e-12)
})

test_that("cpp_splitwise_info_batch: polytomy trees", {
  poly1 <- ape::read.tree(text = "((t1,t2,t3),(t4,t5,t6));")
  poly2 <- ape::read.tree(text = "((t1,t2),(t3,t4,t5,t6));")
  trees <- structure(list(poly1, poly2), class = "multiPhylo")
  splits_list <- as.Splits(trees)
  nTip <- 6L
  batch <- TreeDist:::cpp_splitwise_info_batch(splits_list, nTip)
  r_ref <- vapply(trees, SplitwiseInfo, double(1))
  expect_equal(batch, unname(r_ref), tolerance = 1e-8)
})


# KendallColijn batch paths ----

test_that("KC all-pairs uses pair_diff_euclidean", {
  trees <- ape::as.phylo(0:4, tipLabels = tips8)
  kc <- KendallColijn(trees)
  expect_s3_class(kc, "dist")
  expect_equal(attr(kc, "Size"), 5L)
  m <- as.matrix(kc)
  expect_equal(m[2, 1], KendallColijn(trees[[1]], trees[[2]]),
               tolerance = 1e-10)
  expect_equal(m[4, 3], KendallColijn(trees[[3]], trees[[4]]),
               tolerance = 1e-10)
})

test_that("KC cross-pairs uses vec_diff_euclidean", {
  cross <- KendallColijn(tA8, tB8)
  expect_equal(dim(cross), c(3L, 3L))
  expect_equal(cross[1, 1], KendallColijn(tA8[[1]], tB8[[1]]),
               tolerance = 1e-10)
  expect_equal(cross[2, 3], KendallColijn(tA8[[2]], tB8[[3]]),
               tolerance = 1e-10)
})

test_that("KC all-pairs with SplitVector", {
  trees <- ape::as.phylo(0:3, tipLabels = tips8)
  kc <- KendallColijn(trees, Vector = SplitVector)
  expect_s3_class(kc, "dist")
  m <- as.matrix(kc)
  expect_equal(m[2, 1],
               KendallColijn(trees[[1]], trees[[2]], Vector = SplitVector),
               tolerance = 1e-10)
})

test_that("KC cross-pairs with SplitVector", {
  cross <- KendallColijn(tA8, tB8, Vector = SplitVector)
  expect_equal(dim(cross), c(3L, 3L))
  expect_equal(cross[1, 1],
               KendallColijn(tA8[[1]], tB8[[1]], Vector = SplitVector),
               tolerance = 1e-10)
})


# Large-tree cross-pairs (multi-bin >64 tips) ----

test_that("Cross-pairs with >64-tip trees exercise multi-bin path", {
  skip_on_cran()
  tips100 <- paste0("t", seq_len(100))
  tA100 <- ape::as.phylo(0:2, tipLabels = tips100)
  tB100 <- ape::as.phylo(3:5, tipLabels = tips100)

  # CID cross-pairs
  cid <- ClusteringInfoDistance(tA100, tB100)
  expect_equal(dim(cid), c(3L, 3L))
  expect_equal(cid[1, 1],
               ClusteringInfoDistance(tA100[[1]], tB100[[1]]),
               tolerance = 1e-10)

  # MSD cross-pairs
  msd <- MatchingSplitDistance(tA100, tB100)
  expect_equal(dim(msd), c(3L, 3L))
  expect_equal(msd[1, 1],
               MatchingSplitDistance(tA100[[1]], tB100[[1]]),
               tolerance = 1e-10)

  # DPI cross-pairs
  dpi <- DifferentPhylogeneticInfo(tA100, tB100)
  expect_equal(dim(dpi), c(3L, 3L))
  expect_equal(dpi[1, 1],
               DifferentPhylogeneticInfo(tA100[[1]], tB100[[1]]),
               tolerance = 1e-10)
})

test_that("All-pairs entropy batch with >64-tip trees", {
  skip_on_cran()
  tips100 <- paste0("t", seq_len(100))
  trees100 <- ape::as.phylo(0:2, tipLabels = tips100)
  splits_list <- as.Splits(trees100)

  ce <- TreeDist:::cpp_clustering_entropy_batch(splits_list, 100L)
  ce_ref <- vapply(trees100, ClusteringEntropy, double(1))
  expect_equal(ce, unname(ce_ref), tolerance = 1e-12)

  si <- TreeDist:::cpp_splitwise_info_batch(splits_list, 100L)
  si_ref <- vapply(trees100, SplitwiseInfo, double(1))
  expect_equal(si, unname(si_ref), tolerance = 1e-8)
})


# SPI/MSI all-pairs batch dispatch via .SplitDistanceAllPairs ----

test_that("SPI all-pairs batch agrees with per-pair", {
  batch <- SharedPhylogeneticInfo(trees20)
  m <- as.matrix(batch)
  expect_equal(m[2, 1],
               SharedPhylogeneticInfo(trees20[[1]], trees20[[2]]),
               tolerance = 1e-10)
})

test_that("MSI all-pairs batch agrees with per-pair", {
  batch <- MatchingSplitInfo(trees20)
  m <- as.matrix(batch)
  expect_equal(m[2, 1],
               MatchingSplitInfo(trees20[[1]], trees20[[2]]),
               tolerance = 1e-10)
})


# Nye all-pairs batch dispatch ----

test_that("Nye all-pairs batch agrees with per-pair", {
  batch <- NyeSimilarity(trees20)
  m <- as.matrix(batch)
  expect_equal(m[2, 1],
               NyeSimilarity(trees20[[1]], trees20[[2]]),
               tolerance = 1e-10)
})
