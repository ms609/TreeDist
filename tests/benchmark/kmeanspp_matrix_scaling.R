# Benchmark: KMeansPP.matrix O(n^2) materialized-dist vs O(n * dim) on-the-fly
#
# The old KMeansPP.matrix() built the full n x n Euclidean distance matrix
# (as.matrix(dist(x))) purely to seed the k-means++ loop -- O(n^2) time and
# memory.  The new version computes each centre's distance row on the fly in
# O(n * dim), so it scales to n that the dense matrix cannot fit in RAM.
#
# This script (1) confirms the two give identical clusterings under a fixed
# seed, and (2) reports wall-time and peak memory for the new method at
# n = 6400, 20000, 50000.  The old method is run only at n = 6400; at larger n
# the dense matrix is reported analytically (it is infeasible to allocate, and
# attempting ~20 GB could thrash the machine).
#
# Run against an installed build, e.g.:
#   R CMD INSTALL --library=.dev-kmpp .
#   Rscript -e ".libPaths(c('.dev-kmpp', .libPaths())); source('tests/benchmark/kmeanspp_matrix_scaling.R')"

library(TreeDist)

# ---- old implementation (materialized n x n distance matrix) ----------------

old_kmeanspp_matrix <- function(x, k = 2, nstart = 10, ...) {
  if (k < 2) {
    return(kmeans(x, centers = k, ...))
  }
  n <- dim(x)[[1]]
  ret <- list(tot.withinss = Inf)
  d <- as.matrix(dist(x))
  for (start in seq_len(nstart)) {
    centres <- integer(k)
    centres[1L] <- sample.int(n, 1L)
    min_d <- d[centres[1L], ]
    for (i in 2:k) {
      centres[i] <- sample.int(n, 1L, prob = min_d ^ 2)
      min_d <- pmin.int(min_d, d[centres[i], ])
    }
    proposal <- kmeans(x, centers = x[centres, ], ...)
    if (proposal[["tot.withinss"]] < ret[["tot.withinss"]]) {
      ret <- proposal
    }
  }
  ret
}

# ---- measurement helpers ----------------------------------------------------

# Peak memory via gc(): reset the high-water mark, run, read "max used (Mb)"
# (column 6) summed over the Ncells and Vcells rows.
measure <- function(thunk) {
  invisible(gc(reset = TRUE, full = TRUE))
  elapsed <- system.time(value <- thunk())[["elapsed"]]
  g <- gc(full = TRUE)
  list(value = value, elapsed = elapsed, peakMb = sum(g[, 6]))
}

dense_matrix_gb <- function(n) as.double(n) ^ 2 * 8 / 1e9  # one double n x n matrix

# ---- identity checks --------------------------------------------------------

check_identity <- function(n, dim, k, nstart, seed = 1L, runSeed = 42L) {
  set.seed(seed)
  x <- matrix(rnorm(n * dim), ncol = dim)
  set.seed(runSeed)
  old <- suppressWarnings(old_kmeanspp_matrix(x, k = k, nstart = nstart,
                                              iter.max = 100L))
  set.seed(runSeed)
  new <- suppressWarnings(KMeansPP(x, k = k, nstart = nstart, iter.max = 100L))
  ok <- identical(old[["cluster"]], new[["cluster"]]) &&
    identical(old[["tot.withinss"]], new[["tot.withinss"]])
  cat(sprintf("  n=%5d dim=%2d k=%2d: %s (tot.withinss old=%.6g new=%.6g)\n",
              n, dim, k, if (ok) "IDENTICAL" else "*** DIFFERS ***",
              old[["tot.withinss"]], new[["tot.withinss"]]))
  ok
}

cat("Identity (old materialized-dist vs new on-the-fly), fixed seed:\n")
ok1 <- check_identity(n = 2000, dim =  5, k = 10, nstart = 10)
ok2 <- check_identity(n = 3000, dim = 36, k = 20, nstart =  5)
stopifnot(ok1, ok2)

# ---- timing + peak memory ---------------------------------------------------

K <- 10L
NSTART <- 10L
cat(sprintf("\nTiming + peak memory (k=%d, nstart=%d, dim=36):\n", K, NSTART))
cat(sprintf("  %-7s | %-26s | %-26s | %s\n",
            "n", "old (materialized dist)", "new (on-the-fly)", "speed-up"))

run_size <- function(n, runOld) {
  set.seed(99L)
  x <- matrix(rnorm(n * 36), ncol = 36)

  newRes <- measure(function() {
    set.seed(7L)
    suppressWarnings(KMeansPP(x, k = K, nstart = NSTART, iter.max = 100L))
  })

  if (isTRUE(runOld)) {
    oldRes <- measure(function() {
      set.seed(7L)
      suppressWarnings(old_kmeanspp_matrix(x, k = K, nstart = NSTART,
                                           iter.max = 100L))
    })
    oldStr <- sprintf("%6.2f s / %7.0f MB", oldRes[["elapsed"]],
                      oldRes[["peakMb"]])
    speed <- sprintf("%.1fx", oldRes[["elapsed"]] / newRes[["elapsed"]])
  } else {
    oldStr <- sprintf("~%.1f GB matrix (skipped)", dense_matrix_gb(n))
    speed <- "n/a"
  }
  cat(sprintf("  %-7d | %-26s | %6.2f s / %7.0f MB | %s\n",
              n, oldStr, newRes[["elapsed"]], newRes[["peakMb"]], speed))
}

run_size( 6400L, runOld = TRUE)
run_size(20000L, runOld = FALSE)  # old dense matrix ~3 GB: reported analytically
run_size(50000L, runOld = FALSE)  # old dense matrix ~20 GB: infeasible

cat("\nNote: the old method must allocate an n x n double matrix",
    "(n^2 * 8 bytes),\n   ~3.2 GB at n=20000 and ~20 GB at n=50000",
    "(peak ~1.5x during as.matrix(dist())),\n   so it is effectively unusable",
    "there; the new method's memory grows ~linearly in n.\n")
