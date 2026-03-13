# benchmark/vtune-driver.R
#
# Driver script for VTune hotspot collection.
# Run via:
#   vtune -collect hotspots -knob sampling-mode=hw -result-dir vtune-current \
#         -- Rscript benchmark/vtune-driver.R
#
# The script exercises the two main LAP-heavy distance functions:
#   - ClusteringInfoDistance  (CID): spi_overlap -> CostMatrix -> LAP
#   - PhylogeneticInfoDistance (PID): spi_overlap -> CostMatrix -> LAP
#
# Each workload loops for ~30 s so VTune accumulates enough samples.

vtune_lib <- normalizePath("C:/Users/pjjg18/GitHub/TreeDist/.vtune-lib",
                           mustWork = TRUE)
.libPaths(c(vtune_lib, .libPaths()))

suppressPackageStartupMessages({
  library(TreeTools)
  library(TreeDist)
})

cat("TreeDist version:", as.character(packageVersion("TreeDist")),
    "from:", find.package("TreeDist"), "\n")

set.seed(4271)
tr50  <- as.phylo(0:99,  tipLabels = paste0("t", seq_len(50)))   # 100 trees
tr200 <- as.phylo(0:39,  tipLabels = paste0("t", seq_len(200)))  # 40 trees

TARGET_SECS <- 30

# ---  CID 50-tip  ------------------------------------------------------------
cat("\nWarming up CID 50-tip ...\n")
invisible(ClusteringInfoDistance(tr50))
invisible(ClusteringInfoDistance(tr50))

cat("Running CID 50-tip for ~", TARGET_SECS, "s ...\n")
t0 <- proc.time()[["elapsed"]]
n_iter <- 0L
while (proc.time()[["elapsed"]] - t0 < TARGET_SECS) {
  invisible(ClusteringInfoDistance(tr50))
  n_iter <- n_iter + 1L
}
cat("  iterations:", n_iter, "\n")

# ---  CID 200-tip  -----------------------------------------------------------
cat("\nWarming up CID 200-tip ...\n")
invisible(ClusteringInfoDistance(tr200))

cat("Running CID 200-tip for ~", TARGET_SECS, "s ...\n")
t0 <- proc.time()[["elapsed"]]
n_iter <- 0L
while (proc.time()[["elapsed"]] - t0 < TARGET_SECS) {
  invisible(ClusteringInfoDistance(tr200))
  n_iter <- n_iter + 1L
}
cat("  iterations:", n_iter, "\n")

# ---  PID 50-tip  ------------------------------------------------------------
cat("\nRunning PID 50-tip for ~", TARGET_SECS, "s ...\n")
invisible(PhylogeneticInfoDistance(tr50))
t0 <- proc.time()[["elapsed"]]
n_iter <- 0L
while (proc.time()[["elapsed"]] - t0 < TARGET_SECS) {
  invisible(PhylogeneticInfoDistance(tr50))
  n_iter <- n_iter + 1L
}
cat("  iterations:", n_iter, "\n")

cat("\nDone.\n")
