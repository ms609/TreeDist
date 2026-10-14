# The portable popcount path is what runs on x86-64 CPUs without POPCNT
# (ms609/TreeTools#287). CI machines all have POPCNT, so force that path and
# check that every kernel that uses it returns identical results.

test_that("Portable popcount path matches the POPCNT path", {
  library("TreeTools", quietly = TRUE)
  forceSW <- TreeDist:::cpp_popcount_force_sw
  expect_false(forceSW(FALSE))
  on.exit(forceSW(FALSE), add = TRUE)

  metrics <- list(ClusteringInfoDistance, MutualClusteringInfo,
                  PhylogeneticInfoDistance, SharedPhylogeneticInfo,
                  MatchingSplitDistance, MatchingSplitInfoDistance,
                  InfoRobinsonFoulds, JaccardRobinsonFoulds, NyeSimilarity,
                  TransferDist)

  Scores <- function(nTip) {
    set.seed(nTip)
    trees <- structure(list(BalancedTree(nTip), PectinateTree(nTip),
                            RandomTree(nTip, root = TRUE),
                            RandomTree(nTip, root = TRUE),
                            RandomTree(nTip, root = TRUE)),
                       class = "multiPhylo")
    # Trees that share most splits take the pinned SPI and exact-match paths
    similar <- as.phylo(c(0:3, 40), nTip)
    sp1 <- as.Splits(trees[[1]])
    sp3 <- as.Splits(trees[[3]], trees[[1]])
    c(
      lapply(metrics, function(Metric) Metric(trees)),
      lapply(metrics, function(Metric) Metric(trees[1:2], trees[3:5])),
      lapply(metrics, function(Metric) Metric(trees[[1]], trees[[3]])),
      list(
        JaccardRobinsonFoulds(trees, allowConflict = FALSE),
        JaccardRobinsonFoulds(trees[[1]], trees[[3]], allowConflict = FALSE),
        TreeDist:::mismatch_size(sp1, sp3),
        TreeDist:::confusion(sp1, sp3),
        TreeDist:::cpp_mci_impl_score(sp1, sp3, nTip),
        TreeDist:::cpp_shared_phylo(sp1, sp3, nTip, force_slow = FALSE),
        TreeDist:::cpp_shared_phylo(sp1, sp3, nTip, force_slow = TRUE),
        suppressWarnings(TransferConsensus(trees)),
        SharedPhylogeneticInfo(similar),
        SharedPhylogeneticInfo(similar[1:2], similar[3:5]),
        ClusteringInfoDistance(similar),
        MatchingSplitInfoDistance(similar),
        SPRDist(trees[1:3], method = "Rogue"),
        SPRDist(trees[[1]], trees[[3]]),
        HierarchicalMutualInfo(trees[[1]], trees[[3]])
      )
    )
  }

  for (nTip in c(20, 70, 140)) {
    hw <- Scores(nTip)
    expect_false(forceSW(TRUE))
    sw <- Scores(nTip)
    expect_true(forceSW(FALSE))
    expect_equal(sw, hw, info = paste(nTip, "leaves"))
  }
})
