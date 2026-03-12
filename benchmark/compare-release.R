# benchmark/compare-release.R
#
# Release-quality benchmarking for pairwise distance functions.
#
# devtools::load_all() appends -UNDEBUG -g -O0 to the compiler flags, which
# overrides the -O2 in ~/.R/Makevars and produces an unrepresentative build.
# This script avoids that by installing the package to a private library via
# install.packages(..., type = "source") and running each benchmark in a
# fresh Rscript subprocess (sidestepping Windows DLL lock).
#
# USAGE (from a fresh R session — do NOT load TreeDist first):
#
#   source("benchmark/compare-release.R")
#
#   # Benchmark the current working tree (label defaults to git description)
#   bench_release()
#
#   # Benchmark a specific label (e.g. before applying a patch)
#   bench_release(label = "baseline")
#
#   # Compare two previously saved results
#   bench_compare("baseline", "scratch-reuse")
#
# Results are saved to benchmark/results/<label>.Rds so they persist across
# sessions and can be compared at any time.
# ---------------------------------------------------------------------------

.REPO   <- normalizePath("C:/Users/pjjg18/GitHub/TreeDist")
.RESDIR <- file.path(.REPO, "benchmark", "results")

# ---------------------------------------------------------------------------
# bench_release(label)
# Installs the current working-tree source to a private temp library,
# runs the pairwise benchmarks in a subprocess, and saves results.
# ---------------------------------------------------------------------------
bench_release <- function(label = NULL) {
  if (is.null(label)) {
    sha   <- tryCatch(
      trimws(system2("git", c("-C", .REPO, "describe", "--always", "--dirty"),
                     stdout = TRUE, stderr = FALSE)),
      error = function(e) "unknown"
    )
    label <- sha
  }

  lib_dir <- file.path(tempdir(), paste0("TreeDist_bench_", label))
  dir.create(lib_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(.RESDIR,  recursive = TRUE, showWarnings = FALSE)
  out_rds <- file.path(.RESDIR, paste0(label, ".Rds"))

  # Remove stale .o / .dll files that devtools::load_all() may have left
  # behind with debug flags (-O0).  Without this, make sees up-to-date
  # timestamps and skips recompilation, producing a "release" install that
  # actually contains -O0 objects.
  stale <- Sys.glob(file.path(.REPO, "src", c("*.o", "*.dll", "*.so")))
  if (length(stale)) file.remove(stale)

  message("── Installing '", label, "' to ", lib_dir, " …")
  install.packages(.REPO, lib = lib_dir, repos = NULL, type = "source",
                   INSTALL_opts = "--no-multiarch", quiet = FALSE)

  # Run benchmarks in a subprocess so we get a clean DLL load and no
  # devtools debug flags pollute the timing.
  script <- tempfile(fileext = ".R")
  writeLines(c(
    sprintf('.libPaths(c("%s", .libPaths()))', gsub("\\\\", "/", lib_dir)),
    'suppressPackageStartupMessages({',
    '  library(TreeTools)',
    '  library(TreeDist)',
    '})',
    'set.seed(9137)',
    'trees50  <- as.phylo(0:99,  tipLabels = paste0("t", seq_len(50)))',
    'trees200 <- as.phylo(0:39,  tipLabels = paste0("t", seq_len(200)))',
    'b_cid_50  <- bench::mark(ClusteringInfoDist(trees50),  min_iterations = 5)',
    'b_cid_200 <- bench::mark(ClusteringInfoDist(trees200), min_iterations = 5)',
    'b_msd_50  <- bench::mark(MatchingSplitDistance(trees50),  min_iterations = 5)',
    'b_msd_200 <- bench::mark(MatchingSplitDistance(trees200), min_iterations = 5)',
    sprintf('saveRDS(list(cid_50=b_cid_50, cid_200=b_cid_200, msd_50=b_msd_50, msd_200=b_msd_200), "%s")',
            gsub("\\\\", "/", out_rds))
  ), script)

  message("── Running benchmarks in subprocess …")
  ret <- system2("Rscript", script)
  if (ret != 0L) stop("Benchmark subprocess failed (exit code ", ret, ")")

  message("── Results saved to ", out_rds)
  invisible(readRDS(out_rds))
}

# ---------------------------------------------------------------------------
# bench_compare(label_a, label_b)
# Loads two saved result sets and prints a side-by-side comparison table.
# ---------------------------------------------------------------------------
bench_compare <- function(label_a, label_b) {
  load_res <- function(label) {
    path <- file.path(.RESDIR, paste0(label, ".Rds"))
    if (!file.exists(path))
      stop("No results found for label '", label, "': ", path)
    readRDS(path)
  }

  a <- load_res(label_a)
  b <- load_res(label_b)

  # bench_mark medians are bench_time objects (stored in seconds)
  to_ms <- function(bm) as.numeric(bm$median) * 1e3

  keys <- c("cid_50", "cid_200", "msd_50", "msd_200")
  desc <- c(
    "ClusteringInfoDist  — 100 trees × 50 tips  (4 950 pairs)",
    "ClusteringInfoDist  — 40 trees  × 200 tips  (780 pairs)",
    "MatchingSplitDist   — 100 trees × 50 tips  (4 950 pairs)",
    "MatchingSplitDist   — 40 trees  × 200 tips  (780 pairs)"
  )

  cat(sprintf("\n%-55s  %9s  %9s  %8s\n",
              "Scenario", label_a, label_b, "speedup"))
  cat(strrep("─", 85), "\n")

  for (i in seq_along(keys)) {
    ma <- to_ms(a[[keys[i]]])
    mb <- to_ms(b[[keys[i]]])
    cat(sprintf("%-55s  %7.1f ms  %7.1f ms  %+6.1f%%\n",
                desc[i], ma, mb, (ma - mb) / ma * 100))
  }
  cat("\nNote: positive speedup = ", label_b, " is faster\n", sep = "")
  invisible(list(a = a, b = b))
}
