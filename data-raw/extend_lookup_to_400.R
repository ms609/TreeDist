# Extend TreeDist::randomTreeDistances to cover n=301..400 for CID and PID.
#
# Companion to `extend_lookup_to_300.R`: identical logic with new_ns <- 301:400.
# Uses the *corrected* TreeTools::RandomTree (acid-aphid branch worktree) so
# the new cells are drawn from the uniform distribution over labelled rooted
# binary trees.
#
# ~60-80 min on a typical laptop. SE/mean stays under 0.0002 throughout.

suppressPackageStartupMessages({
  devtools::load_all("C:/Users/pjjg18/GitHub/worktrees/acid-aphid-TreeDist", quiet = TRUE)
  devtools::load_all("C:/Users/pjjg18/GitHub/worktrees/acid-aphid-TreeTools", quiet = TRUE)
})

set.seed(20260514L)

current <- TreeDist::randomTreeDistances
cat("Current rows: n=", range(current$n), "\n", sep = "")

# Refresh 201..300 (previously built under buggy RandomTree) AND extend to 400.
new_ns <- 201:400
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
  flush.console()
  out[[i]] <- data.frame(
    n = n,
    cid_mean = mean(cid), cid_sd = sd(cid), n_repls_cid = nReps,
    pid_mean = mean(pid), pid_sd = sd(pid), n_repls_pid = nReps
  )
}

new_rows <- do.call(rbind, out)
# Replace any existing rows with n in new_ns (e.g. stale n=201..300 from the
# previous extend_to_300 run) before appending.
keep <- !(current$n %in% new_ns)
randomTreeDistances <- rbind(current[keep, , drop = FALSE], new_rows)
randomTreeDistances <- randomTreeDistances[order(randomTreeDistances$n), ]
rownames(randomTreeDistances) <- NULL

pkg_root <- "C:/Users/pjjg18/GitHub/worktrees/acid-aphid-TreeDist"
save_path <- file.path(pkg_root, "data", "randomTreeDistances.rda")
save(randomTreeDistances, file = save_path, compress = "xz")
cat("Saved to:", save_path, "\n")
cat("Final rows: n=", range(randomTreeDistances$n), "  count=", nrow(randomTreeDistances), "\n", sep = "")
cat("File size:", file.info(save_path)$size, "bytes\n")
