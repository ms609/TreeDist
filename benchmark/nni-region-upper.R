# Cost and tightness of NNIDist(regionTips = ...), the region-by-region
# exact upper bound.  Run from the package root against a release install:
#   R CMD INSTALL --library=.agent-lib .
#   Rscript benchmark/nni-region-upper.R
# Not a CI benchmark (no `bench-` prefix).

.libPaths(c(".agent-lib", .libPaths()))
library("TreeTools", quietly = TRUE)
library("TreeDist")
.NNIExact <- TreeDist:::.NNIExact

# Pair of trees on n leaves sharing no splits: a single region of n - 3 edges
.DisjointPair <- function(n) {
  repeat {
    pair <- list(RandomTree(n, root = FALSE), RandomTree(n, root = FALSE))
    pair[[2]] <- RenumberTips(pair[[2]], pair[[1]])
    if (NNIDist(pair[[1]], pair[[2]])[["lower"]] == n - 3) return(pair)
  }
}

.Perturb <- function(tree, moves) {
  out <- phangorn::rNNI(tree, moves = moves)
  out[["edge.length"]] <- NULL
  out
}

.Summarize <- function(pairs, exact = TRUE) {
  res <- t(vapply(pairs, function(pair) {
    d <- NNIDist(pair[[1]], pair[[2]], regionTips = 16L)
    c(d[c("lower", "best_upper", "region_upper")],
      exact = if (exact) .NNIExact(pair[[1]], pair[[2]]) else NA_integer_)
  }, double(4)))
  c(colMeans(res),
    `region < best` = mean(res[, "region_upper"] < res[, "best_upper"]),
    `region == exact` = mean(res[, "region_upper"] == res[, "exact"]))
}

# 1. Cost of one region against its size ####
set.seed(1)
budget <- 2^21 # NNI_EXACT_BUDGET in nni_exact.h
cost <- t(vapply(5:16, function(n) {
  reps <- if (n < 12) 20 else 5
  pairs <- replicate(reps, .DisjointPair(n), simplify = FALSE)
  times <- vapply(pairs, function(pair) {
    system.time(d <- .NNIExact(pair[[1]], pair[[2]], budget = budget))[[3]]
  }, double(1))
  found <- vapply(pairs, function(pair) {
    !is.na(.NNIExact(pair[[1]], pair[[2]], budget = budget))
  }, logical(1))
  c(leaves = n, medianSeconds = median(times), maxSeconds = max(times),
    solved = mean(found))
}, double(4)))
print(cost, digits = 3)

# 2. Tightness at n = 7, 8 against the exact distance ####
set.seed(2)
tight <- do.call(rbind, lapply(7:8, function(n) {
  random <- replicate(400, {
    pair <- list(RandomTree(n, root = FALSE), RandomTree(n, root = FALSE))
    pair[[2]] <- RenumberTips(pair[[2]], pair[[1]])
    pair
  }, simplify = FALSE)
  near <- lapply(1:4, function(moves) replicate(200, {
    tree <- RandomTree(n, root = FALSE)
    list(tree, .Perturb(tree, moves))
  }, simplify = FALSE))
  rbind(
    cbind(n = n, moves = NA, t(.Summarize(random))),
    do.call(rbind, lapply(1:4, function(moves) {
      cbind(n = n, moves = moves, t(.Summarize(near[[moves]])))
    }))
  )
}))
print(tight, digits = 3)

# 3. Larger trees ####
set.seed(3)
large <- do.call(rbind, lapply(c(12, 16, 25, 50, 100), function(n) {
  do.call(rbind, lapply(c(2, 5, 10, 20, NA), function(moves) {
    pairs <- replicate(50, {
      tree <- RandomTree(n, root = FALSE)
      other <- if (is.na(moves)) {
        RenumberTips(RandomTree(n, root = FALSE), tree)
      } else {
        .Perturb(tree, moves)
      }
      list(tree, other)
    }, simplify = FALSE)
    elapsed <- system.time(s <- .Summarize(pairs, exact = FALSE))[[3]]
    cbind(n = n, moves = moves, t(s[1:4]), `region < best` = s[[5]],
          msPerPair = 1000 * elapsed / length(pairs))
  }))
}))
print(large, digits = 3)

# 4. Near-identical pairs where the exact distance is still searchable ####
set.seed(4)
mid <- do.call(rbind, lapply(c(12, 14, 16), function(n) {
  do.call(rbind, lapply(c(2, 4, 6), function(moves) {
    pairs <- replicate(30, {
      tree <- RandomTree(n, root = FALSE)
      list(tree, .Perturb(tree, moves))
    }, simplify = FALSE)
    cbind(n = n, moves = moves, t(.Summarize(pairs)))
  }))
}))
print(mid, digits = 3)

saveRDS(
  list(cost = cost, tight = tight, large = large, mid = mid),
  file.path(Sys.getenv("NNI_BENCH_OUT", tempdir()), "nni-region-upper.Rds")
)
