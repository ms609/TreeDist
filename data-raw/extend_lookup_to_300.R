# Extend TreeDist::randomTreeDistances to cover n=201..300 for CID and PID.
#
# Companion to `build_random_tree_distances.R`: that script reads the
# TreeDistData reference array (n = 4..200); this one tops the slim object
# up to n=300 via direct Monte Carlo (2000 reps per n, CID+PID only) so the
# package can serve adjusted distances for tree pairs with up to 300 leaves
# without falling back to on-the-fly MC.
#
# Rerun whenever the lookup range needs to grow further; ~45 min for n=201..300
# on a typical laptop. SE/mean stays well under 0.0002 throughout this range.

suppressPackageStartupMessages({
  devtools::load_all("C:/Users/pjjg18/GitHub/worktrees/acid-aphid-TreeDist", quiet = TRUE)
  library("TreeTools")
})

set.seed(20260513L)

current <- TreeDist::randomTreeDistances
cat("Current rows: n=", range(current$n), "\n", sep = "")

new_ns <- 201:300
nReps <- 2000L
out <- vector("list", length(new_ns))

for (i in seq_along(new_ns)) {
  n <- new_ns[i]
  tic <- Sys.time()
  cid <- numeric(nReps)
  pid <- numeric(nReps)
  for (r in seq_len(nReps)) {
    t1 <- RandomTree(n, root = TRUE)
    t2 <- RandomTree(n, root = TRUE)
    cid[r] <- ClusteringInfoDistance(t1, t2, normalize = TRUE)
    pid[r] <- PhylogeneticInfoDistance(t1, t2, normalize = TRUE)
  }
  toc <- Sys.time()
  cid_se_rel <- (sd(cid) / sqrt(nReps)) / mean(cid)
  pid_se_rel <- (sd(pid) / sqrt(nReps)) / mean(pid)
  cat(sprintf("n=%3d  wall=%5.1fs  cid_mean=%.4f sd=%.4f SE/mean=%.5f  pid_mean=%.4f sd=%.4f SE/mean=%.5f\n",
              n, as.numeric(toc - tic, units = "secs"),
              mean(cid), sd(cid), cid_se_rel,
              mean(pid), sd(pid), pid_se_rel))
  out[[i]] <- data.frame(
    n = n,
    cid_mean = mean(cid), cid_sd = sd(cid), n_repls_cid = nReps,
    pid_mean = mean(pid), pid_sd = sd(pid), n_repls_pid = nReps
  )
}

new_rows <- do.call(rbind, out)
randomTreeDistances <- rbind(current, new_rows)
randomTreeDistances <- randomTreeDistances[order(randomTreeDistances$n), ]
rownames(randomTreeDistances) <- NULL

pkg_root <- "C:/Users/pjjg18/GitHub/worktrees/acid-aphid-TreeDist"
save_path <- file.path(pkg_root, "data", "randomTreeDistances.rda")
save(randomTreeDistances, file = save_path, compress = "xz")
cat("Saved to:", save_path, "\n")
cat("Final rows: n=", range(randomTreeDistances$n), "  count=", nrow(randomTreeDistances), "\n", sep = "")
cat("File size:", file.info(save_path)$size, "bytes\n")
