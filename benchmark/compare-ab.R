# benchmark/compare-ab.R
#
# A/B benchmarking: install a "reference" copy of the package as TreeDistRef,
# then compare it against the current working-tree source in a single process.
#
# Both builds use release flags (PKG_CXXFLAGS from Makevars.win, plus
# ~/.R/Makevars — NO devtools -O0 override).  Both are loaded into the same
# Rscript subprocess, so CPU state / thermals / background load affect each
# measurement equally.
#
# USAGE:
#
#   source("benchmark/compare-ab.R")
#
#   # 1. On the code you want as baseline:
#   install_ref()
#
#   # 2. Make your changes, then:
#   bench_ab()        # installs dev, loads both, compares in subprocess
#
# ---------------------------------------------------------------------------

.REPO     <- normalizePath("C:/Users/pjjg18/GitHub/TreeDist")
.REF_LIB  <- file.path(.REPO, "benchmark", "ref_lib")
.RESDIR   <- file.path(.REPO, "benchmark", "results")

# ---------------------------------------------------------------------------
# install_ref()
#
# Copies the current source to a temp directory, renames the package to
# TreeDistRef (so both can coexist in the same R session), and installs
# to benchmark/ref_lib/.
# ---------------------------------------------------------------------------
install_ref <- function() {
  pkg <- file.path(tempdir(), "TreeDistRef_src")
  if (dir.exists(pkg)) unlink(pkg, recursive = TRUE)

  # Copy only the directories R CMD INSTALL needs (avoids copying .git/,
  # benchmark/ref_lib/, and other large non-source artifacts).
  dir.create(pkg, recursive = TRUE)
  for (d in c("R", "src", "inst", "man", "vignettes", "data", "data-raw")) {
    src_d <- file.path(.REPO, d)
    if (dir.exists(src_d)) file.copy(src_d, pkg, recursive = TRUE)
  }
  for (f in c("DESCRIPTION", "NAMESPACE", "LICENSE", "LICENSE.md")) {
    src_f <- file.path(.REPO, f)
    if (file.exists(src_f)) file.copy(src_f, pkg)
  }

  # --- Rename package: TreeDist -> TreeDistRef ---

  # DESCRIPTION
  rewrite <- function(path, pattern, replacement) {
    lines <- readLines(path, warn = FALSE)
    lines <- gsub(pattern, replacement, lines)
    writeLines(lines, path)
  }
  rewrite(file.path(pkg, "DESCRIPTION"),
          "^Package: TreeDist$", "Package: TreeDistRef")

  # NAMESPACE  (useDynLib, importFrom, etc.)
  rewrite(file.path(pkg, "NAMESPACE"),
          "\\bTreeDist\\b", "TreeDistRef")

  # R/RcppExports.R  — registered routine symbols: _TreeDist_fn -> _TreeDistRef_fn
  rewrite(file.path(pkg, "R", "RcppExports.R"),
          "_TreeDist_", "_TreeDistRef_")

  # src/RcppExports.cpp — routine symbols + DLL init function
  rewrite(file.path(pkg, "src", "RcppExports.cpp"),
          "_TreeDist_", "_TreeDistRef_")
  rewrite(file.path(pkg, "src", "RcppExports.cpp"),
          "R_init_TreeDist", "R_init_TreeDistRef")

  # Clean any stale object files
  stale <- Sys.glob(file.path(pkg, "src", c("*.o", "*.dll", "*.so")))
  if (length(stale)) file.remove(stale)

  # Install
  dir.create(.REF_LIB, recursive = TRUE, showWarnings = FALSE)
  message("── Installing reference (TreeDistRef) to ", .REF_LIB, " …")
  install.packages(pkg, lib = .REF_LIB, repos = NULL, type = "source",
                   INSTALL_opts = "--no-multiarch", quiet = FALSE)
  message("── Reference installed.")
  invisible(.REF_LIB)
}

# ---------------------------------------------------------------------------
# bench_ab()
#
# Installs the current working-tree source (as TreeDist) to a temp library,
# then runs a subprocess that loads both TreeDistRef and TreeDist and
# compares them with bench::mark().
# ---------------------------------------------------------------------------
bench_ab <- function() {
  # Verify reference exists
  if (!file.exists(file.path(.REF_LIB, "TreeDistRef", "NAMESPACE"))) {
    stop("No reference install found. Run install_ref() first.")
  }

  # Clean stale .o/.dll from dev source
  stale <- Sys.glob(file.path(.REPO, "src", c("*.o", "*.dll", "*.so")))
  if (length(stale)) file.remove(stale)

  # Install dev version to a temp lib
  dev_lib <- file.path(tempdir(), "TreeDist_bench_dev")
  dir.create(dev_lib, recursive = TRUE, showWarnings = FALSE)
  message("── Installing dev (TreeDist) to ", dev_lib, " …")
  install.packages(.REPO, lib = dev_lib, repos = NULL, type = "source",
                   INSTALL_opts = "--no-multiarch", quiet = FALSE)

  dir.create(.RESDIR, recursive = TRUE, showWarnings = FALSE)
  out_rds <- file.path(.RESDIR, "ab.Rds")

  # Build subprocess script
  script <- tempfile(fileext = ".R")
  writeLines(c(
    sprintf('.libPaths(c("%s", "%s", .libPaths()))',
            gsub("\\\\", "/", dev_lib),
            gsub("\\\\", "/", .REF_LIB)),
    'suppressPackageStartupMessages({',
    '  library(TreeTools)',
    '  library(TreeDistRef)',
    '  library(TreeDist)',
    '})',
    '',
    'cat("TreeDistRef loaded from:", find.package("TreeDistRef"), "\\n")',
    'cat("TreeDist    loaded from:", find.package("TreeDist"), "\\n\\n")',
    '',
    'set.seed(9137)',
    'trees50  <- as.phylo(0:99,  tipLabels = paste0("t", seq_len(50)))',
    'trees200 <- as.phylo(0:39,  tipLabels = paste0("t", seq_len(200)))',
    '',
    '# --- Correctness check ---',
    'ref50  <- TreeDistRef::ClusteringInfoDistance(trees50)',
    'dev50  <- TreeDist::ClusteringInfoDistance(trees50)',
    'ref200 <- TreeDistRef::ClusteringInfoDistance(trees200)',
    'dev200 <- TreeDist::ClusteringInfoDistance(trees200)',
    'cat("Max |ref - dev|  50-tip:", max(abs(as.numeric(ref50)  - as.numeric(dev50))),  "\\n")',
    'cat("Max |ref - dev| 200-tip:", max(abs(as.numeric(ref200) - as.numeric(dev200))), "\\n\\n")',
    '',
    '# --- Benchmarks ---',
    'b_cid_50 <- bench::mark(',
    '  ref = TreeDistRef::ClusteringInfoDistance(trees50),',
    '  dev = TreeDist::ClusteringInfoDistance(trees50),',
    '  min_iterations = 5, check = FALSE',
    ')',
    'cat("CID, 100 trees x 50 tips (4 950 pairs)\\n")',
    'print(b_cid_50[, c("expression", "min", "median", "mem_alloc", "n_itr")])',
    'cat("\\n")',
    '',
    'b_cid_200 <- bench::mark(',
    '  ref = TreeDistRef::ClusteringInfoDistance(trees200),',
    '  dev = TreeDist::ClusteringInfoDistance(trees200),',
    '  min_iterations = 5, check = FALSE',
    ')',
    'cat("CID, 40 trees x 200 tips (780 pairs)\\n")',
    'print(b_cid_200[, c("expression", "min", "median", "mem_alloc", "n_itr")])',
    'cat("\\n")',
    '',
    'b_msd_50 <- bench::mark(',
    '  ref = TreeDistRef::MatchingSplitDistance(trees50),',
    '  dev = TreeDist::MatchingSplitDistance(trees50),',
    '  min_iterations = 5, check = FALSE',
    ')',
    'cat("MSD, 100 trees x 50 tips (4 950 pairs)\\n")',
    'print(b_msd_50[, c("expression", "min", "median", "mem_alloc", "n_itr")])',
    'cat("\\n")',
    '',
    'b_msd_200 <- bench::mark(',
    '  ref = TreeDistRef::MatchingSplitDistance(trees200),',
    '  dev = TreeDist::MatchingSplitDistance(trees200),',
    '  min_iterations = 5, check = FALSE',
    ')',
    'cat("MSD, 40 trees x 200 tips (780 pairs)\\n")',
    'print(b_msd_200[, c("expression", "min", "median", "mem_alloc", "n_itr")])',
    'cat("\\n")',
    '',
    sprintf('saveRDS(list(cid_50 = b_cid_50, cid_200 = b_cid_200, msd_50 = b_msd_50, msd_200 = b_msd_200), "%s")',
            gsub("\\\\", "/", out_rds)),
    'cat("Results saved to", normalizePath("', gsub("\\\\", "/", out_rds), '"), "\\n")'
  ), script)

  message("── Running A/B comparison in subprocess …")
  ret <- system2("Rscript", script)
  if (ret != 0L) stop("Benchmark subprocess failed (exit code ", ret, ")")

  invisible(readRDS(out_rds))
}
