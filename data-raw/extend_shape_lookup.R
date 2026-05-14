# extend_shape_lookup.R -----------------------------------------------------
# (1) Rebuild n=4 cells of inst/extdata/shapeExpectedDistances.rds: the
#     previously committed n=4 values were drawn while
#     PhylogeneticInfoDistance was broken at n=4 (upstream bug, fixed in
#     TreeDist commit 18371ac7 "Small-tree comparison fix"); PID is all 1.0
#     in the stale rows.
# (2) Add n=10 to the shipped shape lookup (nSim = 500 per pair; ~2 hours).
#
# n=11 and n=12 deferred to Hamilton: estimated wall time 8 hours and 15+
# hours respectively at usable nSim.  See PLAN.md "Findings to follow up".

suppressPackageStartupMessages({
  devtools::load_all("C:/Users/pjjg18/GitHub/worktrees/acid-aphid-TreeTools",
                     quiet = TRUE)
  devtools::load_all(".", quiet = TRUE)
})

out_file <- file.path("inst", "extdata", "shapeExpectedDistances.rds")
result_list <- readRDS(out_file)
cat("Loaded:", paste(names(result_list), collapse = ","), "\n")

.ShapePairMC <- function(s1, s2, n, nSim) {
  lab_pool <- as.character(seq_len(n))
  cid_v <- numeric(nSim); pid_v <- numeric(nSim)
  for (i in seq_len(nSim)) {
    t1 <- TreeTools::Preorder(
            TreeTools::RootedTreeWithShape(s1, n, sample(lab_pool)))
    t2 <- TreeTools::Preorder(
            TreeTools::RootedTreeWithShape(s2, n, sample(lab_pool)))
    cid_v[i] <- ClusteringInfoDistance(t1, t2, normalize = TRUE)
    pid_v[i] <- PhylogeneticInfoDistance(t1, t2, normalize = TRUE)
  }
  data.frame(
    shape1 = s1, shape2 = s2,
    mean_cid = mean(cid_v), sd_cid = sd(cid_v),
    mean_pid = mean(pid_v), sd_pid = sd(pid_v),
    n_pairs = nSim
  )
}

BuildOneN <- function(n, nSim) {
  W <- as.integer(TreeTools::NRootedShapes(n))
  pairs_total <- W * (W + 1L) / 2L
  cat(sprintf("[n=%d] W=%d shapes, %d pairs, nSim=%d\n", n, W, pairs_total, nSim))
  t0 <- Sys.time()
  rows <- vector("list", pairs_total)
  k <- 0L
  for (s1 in seq_len(W) - 1L) {
    for (s2 in s1:(W - 1L)) {
      k <- k + 1L
      rows[[k]] <- .ShapePairMC(s1, s2, n, nSim)
    }
    if (W > 20 && s1 %% 10L == 0L && s1 > 0L) {
      cat(sprintf("  s1=%d/%d  k=%d/%d  elapsed=%.1f s\n",
                  s1, W - 1L, k, pairs_total,
                  as.numeric(Sys.time() - t0, units = "secs")))
      flush.console()
    }
  }
  df <- do.call(rbind, rows)
  cat(sprintf("  done in %.1f s\n",
              as.numeric(Sys.time() - t0, units = "secs")))
  df
}

# (1) Rebuild n=4 against fixed PID.
result_list[["4"]] <- BuildOneN(4L, 1000L)
saveRDS(result_list, out_file, compress = TRUE)
cat("After n=4 refresh, file size:", file.info(out_file)$size, "bytes\n")
cat("n=4 row check:\n"); print(result_list[["4"]])

# (2) Add n=10.
if (!"10" %in% names(result_list)) {
  result_list[["10"]] <- BuildOneN(10L, 500L)
  result_list <- result_list[order(as.integer(names(result_list)))]
  saveRDS(result_list, out_file, compress = TRUE)
}

cat("\nFinal keys:", paste(names(result_list), collapse = ","), "\n")
cat("Final file size:", file.info(out_file)$size, "bytes\n")
