# shape_lookup_n.R - compute the shape-conditioned mean CID & PID for every
# pair of rooted shapes at given n.  Args: n nSim out_path
suppressPackageStartupMessages({
  .libPaths(c(Sys.getenv("R_LIBS_USER", "/nobackup/pjjg18/acid-aphid/lib"),
              .libPaths()))
  library("TreeTools")
  library("TreeDist")
})

args <- commandArgs(trailingOnly = TRUE)
n    <- as.integer(args[1])
nSim <- as.integer(args[2])
out_path <- args[3]
stopifnot(!is.na(n), !is.na(nSim), nzchar(out_path))

cat(sprintf("[%s] n=%d nSim=%d out=%s\n",
            format(Sys.time(), "%Y-%m-%d %H:%M:%S"), n, nSim, out_path))

W <- as.integer(NRootedShapes(n))
pairs_total <- W * (W + 1L) / 2L
cat(sprintf("W=%d  pairs=%d\n", W, pairs_total))

ShapePairMC <- function(s1, s2, n, nSim) {
  lab_pool <- as.character(seq_len(n))
  cid_v <- numeric(nSim); pid_v <- numeric(nSim)
  for (i in seq_len(nSim)) {
    t1 <- Preorder(RootedTreeWithShape(s1, n, sample(lab_pool)))
    t2 <- Preorder(RootedTreeWithShape(s2, n, sample(lab_pool)))
    cid_v[i] <- ClusteringInfoDistance(t1, t2, normalize = TRUE)
    pid_v[i] <- PhylogeneticInfoDistance(t1, t2, normalize = TRUE)
  }
  c(mean_cid = mean(cid_v), sd_cid = sd(cid_v),
    mean_pid = mean(pid_v), sd_pid = sd(pid_v))
}

t0 <- Sys.time()
rows <- vector("list", pairs_total)
k <- 0L
for (s1 in seq_len(W) - 1L) {
  for (s2 in s1:(W - 1L)) {
    k <- k + 1L
    stats <- ShapePairMC(s1, s2, n, nSim)
    rows[[k]] <- data.frame(shape1 = s1, shape2 = s2,
                            mean_cid = stats[["mean_cid"]],
                            sd_cid   = stats[["sd_cid"]],
                            mean_pid = stats[["mean_pid"]],
                            sd_pid   = stats[["sd_pid"]],
                            n_pairs  = nSim)
  }
  if (s1 %% 10L == 0L && s1 > 0L) {
    elapsed <- as.numeric(Sys.time() - t0, units = "secs")
    eta <- elapsed * pairs_total / k - elapsed
    cat(sprintf("[%s] s1=%d/%d  k=%d/%d  elapsed=%.1fs  eta=%.1fs\n",
                format(Sys.time(), "%H:%M:%S"),
                s1, W - 1L, k, pairs_total, elapsed, eta))
    flush.console()
  }
}
df <- do.call(rbind, rows)
saveRDS(df, out_path, compress = TRUE)
cat(sprintf("Wrote %s (%.1f KB)\n", out_path, file.info(out_path)$size / 1024))
cat(sprintf("Total wall: %.1f s\n", as.numeric(Sys.time() - t0, units = "secs")))
