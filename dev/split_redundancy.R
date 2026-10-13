# Is a global unique-split index (issue #57) worth building?
#
# For tree sets of varying redundancy, report
#  - split redundancy: unique splits / total splits, and the mean fraction of
#    splits matched exactly between a pair of trees;
#  - wall time of the R-side split conversion versus the C++ all-pairs call;
#  - (optionally) the per-phase breakdown of the C++ call, from a build
#    instrumented to expose `cpp_phase_times()`.
#
# A unique-split index could at best remove the exact-match scan and replace
# the cost-matrix fill with table look-ups; it cannot touch the LAP.  So the
# fill + match share of run time bounds the achievable gain.
#
# Usage: Rscript dev/split_redundancy.R [libPath]
#   where libPath optionally holds a TreeDist build with
#   dev/split_redundancy_phase_timers.patch applied (single-threaded only).

args <- commandArgs(trailingOnly = TRUE)
if (length(args)) .libPaths(c(args[[1]], .libPaths()))
suppressPackageStartupMessages({
  library("TreeTools")
  library("TreeDist")
})
phaseTimes <- if (exists("cpp_phase_times", asNamespace("TreeDist"))) {
  get("cpp_phase_times", asNamespace("TreeDist"))
}
set.seed(5757)

RandomNNIWalk <- function(nTip, nTree, movesPerStep) {
  # Successive trees differ by `movesPerStep` NNI rearrangements:
  # a stand-in for a thinned MCMC sample.
  tr <- RandomTree(nTip, root = FALSE)
  trees <- vector("list", nTree)
  for (i in seq_len(nTree)) {
    tr <- phangorn::rNNI(tr, moves = movesPerStep)
    trees[[i]] <- tr
  }
  structure(trees, class = "multiPhylo")
}

PerturbedSet <- function(nTip, nTree, moves) {
  # Independent perturbations of one base tree: a tight "posterior".
  base <- RandomTree(nTip, root = FALSE)
  structure(lapply(seq_len(nTree), function(i) {
    phangorn::rNNI(base, moves = moves)
  }), class = "multiPhylo")
}

Redundancy <- function(trees, labs) {
  nTip <- length(labs)
  canon <- lapply(as.Splits(trees, tipLabels = labs), function(s) {
    # Canonicalize so tip 1 is always in the split
    apply(as.logical(s), 1, function(x) {
      if (!x[1]) x <- !x
      paste(as.integer(x), collapse = "")
    })
  })
  total <- sum(lengths(canon))
  unique <- length(unique(unlist(canon)))
  pairs <- combn(length(canon), 2)
  if (ncol(pairs) > 2000) pairs <- pairs[, sample.int(ncol(pairs), 2000)]
  shared <- mean(apply(pairs, 2, function(p) {
    length(intersect(canon[[p[1]]], canon[[p[2]]])) / (nTip - 3)
  }))
  c(splits = total, unique = unique, uniqueRatio = unique / total,
    pairShared = shared)
}

Scenario <- function(label, trees) {
  nTip <- NTip(trees[[1]])
  labs <- TipLabels(trees[[1]])
  tConv <- system.time(
    for (i in 1:5) splits <- as.Splits(trees, tipLabels = labs, asSplits = FALSE)
  )[["elapsed"]] / 5
  red <- Redundancy(trees, labs)

  Time <- function(Fn) {
    if (!is.null(phaseTimes)) phaseTimes(TRUE)
    reps <- 2
    t <- system.time(for (i in seq_len(reps)) Fn(splits, nTip, 1L))[["elapsed"]]
    ph <- if (is.null(phaseTimes)) rep(NA, 3) else phaseTimes(TRUE)[1:3] / 1e9
    c(total = t / reps, ph / reps)
  }
  mci <- Time(TreeDist:::cpp_mutual_clustering_all_pairs)
  pid <- Time(TreeDist:::cpp_shared_phylo_all_pairs)
  data.frame(
    scenario = label, nTip = nTip, nTree = length(trees),
    uniqueRatio = round(red[["uniqueRatio"]], 3),
    pairShared = round(red[["pairShared"]], 3),
    convMs = round(1000 * tConv, 1),
    mciMs = round(1000 * mci[["total"]], 1),
    mciMatch = round(mci[[2]] / mci[["total"]], 3),
    mciFill = round(mci[[3]] / mci[["total"]], 3),
    mciLap = round(mci[[4]] / mci[["total"]], 3),
    pidMs = round(1000 * pid[["total"]], 1),
    pidFill = round(pid[[3]] / pid[["total"]], 3),
    pidLap = round(pid[[4]] / pid[["total"]], 3)
  )
}

scenarios <- list(
  list("random", function(n, N) {
    structure(lapply(seq_len(N), function(i) RandomTree(n, root = FALSE)),
              class = "multiPhylo")
  }),
  list("walk-10NNI", function(n, N) RandomNNIWalk(n, N, 10)),
  list("walk-2NNI", function(n, N) RandomNNIWalk(n, N, 2)),
  list("perturb-5NNI", function(n, N) PerturbedSet(n, N, 5))
)
sizes <- list(c(50, 200), c(200, 100), c(1000, 30))

results <- do.call(rbind, lapply(sizes, function(sz) {
  do.call(rbind, lapply(scenarios, function(sc) {
    Scenario(sc[[1]], sc[[2]](sz[[1]], sz[[2]]))
  }))
}))
options(width = 200)
print(results, row.names = FALSE)
