# benchmark/vtune-lap-driver.R
#
# VTune driver for LAPJV standalone + tree distance profiling.
# Run via:
#   "C:/Program Files (x86)/Intel/oneAPI/vtune/2025.10/bin64/vtune.exe" \
#     -collect hotspots -knob sampling-mode=hw -result-dir vtune-lap \
#     -- Rscript benchmark/vtune-lap-driver.R

vtune_lib <- normalizePath("C:/Users/pjjg18/GitHub/TreeDist/.vtune-lib",
                           mustWork = TRUE)
.libPaths(c(vtune_lib, .libPaths()))

suppressPackageStartupMessages({
  library(TreeTools)
  library(TreeDist)
})

cat("TreeDist version:", as.character(packageVersion("TreeDist")),
    "from:", find.package("TreeDist"), "\n")

TARGET_SECS <- 30

# --- LAPJV 1999x1999 (primary target — regression observed here) ---------
set.seed(7391)
mat2000 <- matrix(runif(1999 * 1999), 1999, 1999)

cat("\nWarming up LAPJV 1999x1999 ...\n")
invisible(LAPJV(mat2000))

cat("Running LAPJV 1999x1999 for ~", TARGET_SECS, "s ...\n")
t0 <- proc.time()[["elapsed"]]
n_iter <- 0L
while (proc.time()[["elapsed"]] - t0 < TARGET_SECS) {
  invisible(LAPJV(mat2000))
  n_iter <- n_iter + 1L
}
cat("  iterations:", n_iter, "\n")

# --- LAPJV 400x400 --------------------------------------------------------
set.seed(7392)
mat400 <- matrix(runif(400 * 400), 400, 400)

cat("\nWarming up LAPJV 400x400 ...\n")
invisible(LAPJV(mat400))

cat("Running LAPJV 400x400 for ~", TARGET_SECS, "s ...\n")
t0 <- proc.time()[["elapsed"]]
n_iter <- 0L
while (proc.time()[["elapsed"]] - t0 < TARGET_SECS) {
  invisible(LAPJV(mat400))
  n_iter <- n_iter + 1L
}
cat("  iterations:", n_iter, "\n")

# --- CID on random 200-tip trees (LAP-heavy) ------------------------------
set.seed(7393)
rand200 <- lapply(1:40, function(i) rtree(200))
class(rand200) <- "multiPhylo"

cat("\nWarming up CID random 200-tip ...\n")
invisible(ClusteringInfoDistance(rand200))

cat("Running CID random 200-tip for ~", TARGET_SECS, "s ...\n")
t0 <- proc.time()[["elapsed"]]
n_iter <- 0L
while (proc.time()[["elapsed"]] - t0 < TARGET_SECS) {
  invisible(ClusteringInfoDistance(rand200))
  n_iter <- n_iter + 1L
}
cat("  iterations:", n_iter, "\n")

# --- MSD on random 200-tip trees (LAP-heavy, no exact-match shortcut) -----
cat("\nRunning MSD random 200-tip for ~", TARGET_SECS, "s ...\n")
invisible(MatchingSplitDistance(rand200))
t0 <- proc.time()[["elapsed"]]
n_iter <- 0L
while (proc.time()[["elapsed"]] - t0 < TARGET_SECS) {
  invisible(MatchingSplitDistance(rand200))
  n_iter <- n_iter + 1L
}
cat("  iterations:", n_iter, "\n")

cat("\nDone.\n")
