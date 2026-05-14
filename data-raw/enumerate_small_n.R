# Exact enumeration of E[CID] and E[PID] for small n (4..7).
#
# Background: `randomTreeDistances` is built by Monte Carlo, drawing two
# independent rooted binary trees with `TreeTools::RandomTree(n, root = TRUE)`
# per replicate.  For small n, the number of distinct trees is small enough
# that we can compute the expectation exactly.
#
# Convention: the MC protocol draws two trees *independently*, so a pair can
# coincide (with probability 1/N for N rooted binary trees on n tips).  The
# unbiased exact match is therefore the mean over all `N^2` ordered pairs,
# *including* self-pairs (which contribute zero distance).  Equivalently:
#
#   E_inc = (2 * sum_{i<j} d(t_i, t_j)) / N^2
#
# At n=4 this matters (N=15; include-self mean differs from exclude-self mean
# by a factor of 14/15).  At n=7 (N=10395) the difference is negligible.
#
# Run sequence: this script first (slow, independent), then
# `build_random_tree_distances.R`, which splices these exact values into the
# slim shipped object.  Reads nothing from the package data; writes
# `data-raw/exact_small_n.rds`.
#
# Run interactively from the package root: `source('data-raw/enumerate_small_n.R')`.

suppressPackageStartupMessages({
  library("TreeTools")
})

# Use this package's distance implementations.
devtools::load_all(quiet = TRUE)

enumerate_one <- function(n) {
  nTrees <- TreeTools::NRooted(n)
  trees <- ape::as.phylo(seq_len(nTrees) - 1L, nTip = n)
  message(sprintf("n = %d: %d trees, %.3g unordered off-diagonal pairs",
                  n, nTrees, nTrees * (nTrees - 1) / 2))

  t0 <- Sys.time()
  d_cid <- ClusteringInfoDistance(trees, normalize = TRUE)
  message(sprintf("  CID computed in %.1f s",
                  as.numeric(Sys.time() - t0, units = "secs")))

  t0 <- Sys.time()
  d_pid <- PhylogeneticInfoDistance(trees, normalize = TRUE)
  message(sprintf("  PID computed in %.1f s",
                  as.numeric(Sys.time() - t0, units = "secs")))

  # Include-self mean & sd.  With N^2 ordered pairs, sum_total = 2 * sum(d)
  # (off-diagonals counted twice, diagonals = 0).  Variance is computed about
  # the include-self mean.
  N2 <- nTrees * nTrees

  inc_mean <- function(d) (2 * sum(d)) / N2
  inc_sd <- function(d) {
    m <- inc_mean(d)
    # sum of squares over all N^2 ordered pairs = 2 * sum(d^2);
    # diagonals contribute 0.  Population variance (denom N^2) is the natural
    # exact-enumeration analogue of the MC `sd` (which uses sample sd over
    # independent draws).  For the small-n table we use the population sd; the
    # discrepancy with the sample-sd MC convention is `sqrt(N^2 / (N^2 - 1))`,
    # i.e. negligible for n >= 5 (N^2 >= 11025).
    var_pop <- (2 * sum(d * d)) / N2 - m * m
    sqrt(pmax(var_pop, 0))
  }

  list(
    n = n,
    nTrees = nTrees,
    cid_mean = inc_mean(d_cid),
    cid_sd = inc_sd(d_cid),
    pid_mean = inc_mean(d_pid),
    pid_sd = inc_sd(d_pid)
  )
}

ns <- 4:7
out <- lapply(ns, enumerate_one)
names(out) <- as.character(ns)

# Flat data.frame for spliceability.
exact_small_n <- do.call(rbind, lapply(out, function(x) {
  data.frame(n = x$n,
             cid_mean = x$cid_mean, cid_sd = x$cid_sd,
             pid_mean = x$pid_mean, pid_sd = x$pid_sd,
             nTrees = x$nTrees,
             stringsAsFactors = FALSE)
}))
row.names(exact_small_n) <- NULL

message("Exact small-n results:")
print(exact_small_n, digits = 8)

saveRDS(exact_small_n, file = "data-raw/exact_small_n.rds")
message("Wrote data-raw/exact_small_n.rds")
