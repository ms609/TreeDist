# Benchmark: streaming pmin vs apply-rebuild in KMeansPP initialisation
#
# Compares old (apply-rebuild) vs new (streaming pmin) seeding loop.
# Also verifies bit-identical output under fixed set.seed().

library(TreeDist)

# ---- helpers ----------------------------------------------------------------

old_seed_loop <- function(d, n, k, nstart) {
  best_centres <- NULL
  best_score <- Inf
  for (s in seq_len(nstart)) {
    centres <- c(sample.int(n, 1L), numeric(k - 1L))
    for (i in 2:k) {
      p <- apply(d[centres, , drop = FALSE], 2, min) ^ 2
      centres[i] <- sample.int(n, 1L, prob = p)
    }
    score <- sum(apply(d[centres, , drop = FALSE], 2, min) ^ 2)
    if (score < best_score) { best_score <- score; best_centres <- centres }
  }
  best_centres
}

new_seed_loop <- function(d, n, k, nstart) {
  best_centres <- NULL
  best_score <- Inf
  for (s in seq_len(nstart)) {
    centres <- integer(k)
    centres[1L] <- sample.int(n, 1L)
    min_d <- d[centres[1L], ]
    for (i in 2:k) {
      p <- min_d ^ 2
      centres[i] <- sample.int(n, 1L, prob = p)
      min_d <- pmin.int(min_d, d[centres[i], ])
    }
    score <- sum(min_d ^ 2)
    if (score < best_score) { best_score <- score; best_centres <- centres }
  }
  best_centres
}

# ---- identity check ---------------------------------------------------------

check_identity <- function(n, k, nstart, seed = 42L) {
  set.seed(seed); x <- matrix(rnorm(n * 5), n, 5)
  d <- as.matrix(dist(x))
  # old loop stores centres as double (c(sample.int, numeric(...))); new as
  # integer. Coerce both to integer before comparing.
  set.seed(seed + 1L); old <- as.integer(old_seed_loop(d, n, k, nstart))
  set.seed(seed + 1L); new <- as.integer(new_seed_loop(d, n, k, nstart))
  identical(old, new)
}

cat("Identity checks:\n")
for (cfg in list(c(200, 10), c(500, 20))) {
  n <- cfg[1]; k <- cfg[2]
  ok <- check_identity(n, k, nstart = 10L)
  cat(sprintf("  N=%d k=%d: %s\n", n, k, if (ok) "PASS" else "FAIL"))
}

# ---- benchmark --------------------------------------------------------------

run_bench <- function(n, k, nstart = 10L, times = 5L) {
  set.seed(99L); x <- matrix(rnorm(n * 5), n, 5)
  d <- as.matrix(dist(x))

  if (requireNamespace("microbenchmark", quietly = TRUE)) {
    mb <- microbenchmark::microbenchmark(
      old = { set.seed(1L); old_seed_loop(d, n, k, nstart) },
      new = { set.seed(1L); new_seed_loop(d, n, k, nstart) },
      times = times
    )
    med <- tapply(mb$time, mb$expr, median)
    ratio <- med[["old"]] / med[["new"]]
    cat(sprintf("  N=%d k=%d: old %.2f ms | new %.2f ms | ratio %.2fx\n",
                n, k,
                med[["old"]] / 1e6,
                med[["new"]] / 1e6,
                ratio))
  } else {
    t_old <- median(replicate(times, {
      set.seed(1L); system.time(old_seed_loop(d, n, k, nstart))["elapsed"]
    }))
    t_new <- median(replicate(times, {
      set.seed(1L); system.time(new_seed_loop(d, n, k, nstart))["elapsed"]
    }))
    ratio <- t_old / t_new
    cat(sprintf("  N=%d k=%d: old %.4f s | new %.4f s | ratio %.2fx\n",
                n, k, t_old, t_new, ratio))
  }
}

cat("\nBenchmark (seeding loop only, nstart=10):\n")
for (cfg in list(c(200, 10), c(500, 20))) {
  run_bench(n = cfg[1], k = cfg[2])
}
