.OracleTrees <- function(nTip) {
  set.seed(nTip)
  trees <- c(
    unclass(ape::as.phylo(sample.int(1e5, 6), nTip)),
    list(
      TreeTools::CollapseNode(ape::as.phylo(1, nTip), nTip + 2:4),
      TreeTools::StarTree(nTip)
    )
  )
  trees[[2]] <- trees[[1]] # Duplicate, so self-distance must be zero
  names(trees) <- paste0("t", seq_along(trees))
  trees
}

.ExpectOracle <- function(trees, Func, ...) {
  oracle <- DistanceOracle(trees, Func, ...)
  for (i in seq_along(trees)) {
    expect_equal(oracle(i), Func(trees[[i]], trees, ...), tolerance = 1e-12)
  }
  invisible(oracle)
}

test_that("DistanceOracle() matches Func(trees[[i]], trees)", {
  skip_on_cran()
  Max <- function(x, y) max(x, y)
  cases <- list(
    list(ClusteringInfoDistance),
    list(InfoRobinsonFoulds),
    list(InfoRobinsonFoulds, similarity = TRUE),
    list(RobinsonFoulds),
    list(RobinsonFoulds, similarity = TRUE),
    list(MatchingSplitDistance),
    list(NyeSimilarity),
    list(NyeSimilarity, normalizeMax = TRUE),
    list(NyeSimilarity, similarity = FALSE),
    list(JaccardRobinsonFoulds),
    list(JaccardRobinsonFoulds, k = 2, allowConflict = FALSE, similarity = TRUE)
  )
  for (nTip in c(9, 140, 330)) {
    trees <- .OracleTrees(nTip)
    for (case in cases) {
      for (normalize in list(FALSE, TRUE, 17, Max)) {
        if (identical(case[[1]], MatchingSplitDistance) &&
            !is.numeric(normalize) && !isFALSE(normalize)) {
          next
        }
        do.call(.ExpectOracle, c(list(trees), case, normalize = normalize))
      }
    }
  }
})

test_that("DistanceOracle() handles Splits, vectors and threads", {
  trees <- .OracleTrees(70)
  splits <- lapply(trees, TreeTools::as.Splits, TreeTools::TipLabels(trees[[1]]))
  oracle <- .ExpectOracle(splits, ClusteringInfoDistance)

  expect_equal(oracle(c(1, 4)),
               cbind(oracle(1), oracle(4)), ignore_attr = TRUE)
  expect_equal(dim(oracle(1:3)), c(length(trees), 3))
  expect_equal(unname(oracle(1)[1:2]), c(0, 0))
  expect_error(oracle(0), "`i` must index trees 1 to 8")
  expect_error(oracle(NA), "`i` must index")
  expect_error(oracle("1"), "`i` must index")
  expect_error(oracle(integer(0)), "`i` must index")

  withr::local_options(mc.cores = 2L)
  .ExpectOracle(trees, MatchingSplitDistance)
})

test_that("DistanceOracle() survives serialization", {
  trees <- .OracleTrees(20)
  oracle <- DistanceOracle(trees, InfoRobinsonFoulds, normalize = TRUE)
  expected <- oracle(3)
  file <- withr::local_tempfile(fileext = ".rds")
  saveRDS(oracle, file)
  restored <- readRDS(file)
  expect_equal(restored(3), expected)
  expect_equal(restored(5), oracle(5))
})

test_that("DistanceOracle() falls back to Func", {
  trees <- .OracleTrees(12)
  bigOracle <- DistanceOracle(trees, RobinsonFoulds, normalize = TRUE,
                              maxBytes = 0)
  expect_equal(bigOracle(3), RobinsonFoulds(trees[[3]], trees, normalize = TRUE))
  expect_equal(dim(bigOracle(1:2)), c(length(trees), 2))

  permuted <- c(trees[1:3], list(TreeTools::BalancedTree(paste0("t", 12:1))))
  .ExpectOracle(permuted, ClusteringInfoDistance)
  mixed <- c(trees[1:3], list(TreeTools::BalancedTree(paste0("u", 1:12))))
  .ExpectOracle(mixed, ClusteringInfoDistance)
  .ExpectOracle(mixed, JaccardRobinsonFoulds, k = 2)

  tiny <- list(TreeTools::BalancedTree(3), TreeTools::PectinateTree(3))
  .ExpectOracle(tiny, RobinsonFoulds)
})

test_that("DistanceOracle() rejects invalid input", {
  trees <- .OracleTrees(8)
  expect_error(DistanceOracle(trees[[1]]), "list of trees")
  expect_error(DistanceOracle(trees, PathDist), "`Func` must be one of")
  expect_error(DistanceOracle(trees, RobinsonFoulds, reportMatching = TRUE),
               "`reportMatching` is not supported")
  expect_error(DistanceOracle(trees, RobinsonFoulds, unknown = 1),
               "unused argument")
  msg <- "Please specify a function to generate a normalizing constant"
  expect_error(MatchingSplitDistance(trees[[1]], trees, normalize = TRUE), msg)
  expect_error(DistanceOracle(trees, MatchingSplitDistance, normalize = TRUE),
               msg)
  expect_error(DistanceOracle(trees, MatchingSplitDistance, normalize = max),
               msg)
})

test_that("cpp_oracle_*() validate input", {
  trees <- .OracleTrees(8)
  splits <- TreeTools::as.Splits(trees, asSplits = FALSE)
  expect_null(cpp_oracle_column(NULL, 1L))
  ptr <- cpp_oracle_prepare(splits, 8L, 1L, 1L)
  expect_gt(attr(ptr, "bytes"), 0)
  expect_error(cpp_oracle_column(ptr, 0L), "`i` must lie between 1 and 8")
  splits[[2]] <- TreeTools::as.Splits(
    TreeTools::BalancedTree(80), asSplits = FALSE
  )
  expect_error(cpp_oracle_prepare(splits, 8L, 1L, 1L), "same number of leaves")
})
