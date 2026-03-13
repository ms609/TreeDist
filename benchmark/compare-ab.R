# benchmark/compare-ab.R
#
# A/B benchmarking: install renamed copies of the package as TreeDistRef and
# TreeDistDev, then compare them in a single Rscript subprocess.
#
# Both builds use release flags (PKG_CXXFLAGS from Makevars.win, plus
# ~/.R/Makevars — NO devtools -O0 override).  Both packages coexist in the
# same process (different DLL names), so CPU state / thermals / background
# load affect each measurement equally.
#
# The rename trick also avoids Windows DLL locks: even if TreeDist is loaded
# in another R session, TreeDistRef and TreeDistDev have different DLL names
# and install without interference.
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
.DEV_LIB  <- file.path(.REPO, "benchmark", "dev_lib")
.RESDIR   <- file.path(.REPO, "benchmark", "results")

# ---------------------------------------------------------------------------
# .install_renamed()
#
# Copies the current source tree, renames the package (DESCRIPTION, NAMESPACE,
# RcppExports symbols, R_init_ entry point), cleans stale object files, and
# installs to `lib_dir`.
# ---------------------------------------------------------------------------
.install_renamed <- function(new_name, lib_dir) {
  pkg <- file.path(tempdir(), paste0(new_name, "_src"))
  if (dir.exists(pkg)) unlink(pkg, recursive = TRUE)

  # Copy only the directories R CMD INSTALL needs (avoids .git/, benchmark/, etc.)
  dir.create(pkg, recursive = TRUE)
  for (d in c("R", "src", "inst", "man", "vignettes", "data", "data-raw")) {
    src_d <- file.path(.REPO, d)
    if (dir.exists(src_d)) file.copy(src_d, pkg, recursive = TRUE)
  }
  for (f in c("DESCRIPTION", "NAMESPACE", "LICENSE", "LICENSE.md")) {
    src_f <- file.path(.REPO, f)
    if (file.exists(src_f)) file.copy(src_f, pkg)
  }

  # --- Rename package: TreeDist -> new_name ---
  rewrite <- function(path, pattern, replacement) {
    lines <- readLines(path, warn = FALSE)
    lines <- gsub(pattern, replacement, lines)
    writeLines(lines, path)
  }

  rewrite(file.path(pkg, "DESCRIPTION"),
          "^Package: TreeDist$", paste0("Package: ", new_name))

  rewrite(file.path(pkg, "NAMESPACE"),
          "\\bTreeDist\\b", new_name)

  rewrite(file.path(pkg, "R", "RcppExports.R"),
          "_TreeDist_", paste0("_", new_name, "_"))

  rewrite(file.path(pkg, "src", "RcppExports.cpp"),
          "_TreeDist_", paste0("_", new_name, "_"))
  rewrite(file.path(pkg, "src", "RcppExports.cpp"),
          "R_init_TreeDist", paste0("R_init_", new_name))

  # Clean stale object files
  stale <- Sys.glob(file.path(pkg, "src", c("*.o", "*.dll", "*.so")))
  if (length(stale)) file.remove(stale)

  # Install
  dir.create(lib_dir, recursive = TRUE, showWarnings = FALSE)
  message("── Installing ", new_name, " to ", lib_dir, " …")
  install.packages(pkg, lib = lib_dir, repos = NULL, type = "source",
                   INSTALL_opts = "--no-multiarch", quiet = FALSE)
  message("── ", new_name, " installed.")
  invisible(lib_dir)
}

# ---------------------------------------------------------------------------
# install_ref() / install_dev()
#
# Install the current working-tree source as TreeDistRef / TreeDistDev.
# Call install_ref() on baseline code, then make changes and call install_dev()
# (or bench_ab(), which calls install_dev() automatically).
# ---------------------------------------------------------------------------
install_ref <- function() .install_renamed("TreeDistRef", .REF_LIB)
install_dev <- function() .install_renamed("TreeDistDev", .DEV_LIB)

# ---------------------------------------------------------------------------
# bench_ab()
#
# Installs the current working-tree source as TreeDistDev, then runs a
# subprocess that loads both TreeDistRef and TreeDistDev and compares
# them with bench::mark().
# ---------------------------------------------------------------------------
bench_ab <- function() {
  if (!file.exists(file.path(.REF_LIB, "TreeDistRef", "NAMESPACE"))) {
    stop("No reference install found. Run install_ref() first.")
  }

  install_dev()

  dir.create(.RESDIR, recursive = TRUE, showWarnings = FALSE)
  out_rds <- file.path(.RESDIR, "ab.Rds")

  # Build subprocess script
  script <- tempfile(fileext = ".R")
  writeLines(c(
    sprintf('.libPaths(c("%s", "%s", .libPaths()))',
            gsub("\\\\", "/", .DEV_LIB),
            gsub("\\\\", "/", .REF_LIB)),
    'suppressPackageStartupMessages({',
    '  library(TreeTools)',
    '  library(TreeDistRef)',
    '  library(TreeDistDev)',
    '})',
    '',
    'cat("TreeDistRef loaded from:", find.package("TreeDistRef"), "\\n")',
    'cat("TreeDistDev loaded from:", find.package("TreeDistDev"), "\\n\\n")',
    '',
    'set.seed(9137)',
    'trees50  <- as.phylo(0:99,  tipLabels = paste0("t", seq_len(50)))',
    'trees200 <- as.phylo(0:39,  tipLabels = paste0("t", seq_len(200)))',
    '',
    '# --- Correctness check ---',
    'ref50  <- TreeDistRef::ClusteringInfoDistance(trees50)',
    'dev50  <- TreeDistDev::ClusteringInfoDistance(trees50)',
    'ref200 <- TreeDistRef::ClusteringInfoDistance(trees200)',
    'dev200 <- TreeDistDev::ClusteringInfoDistance(trees200)',
    'cat("Max |ref - dev|  50-tip:", max(abs(as.numeric(ref50)  - as.numeric(dev50))),  "\\n")',
    'cat("Max |ref - dev| 200-tip:", max(abs(as.numeric(ref200) - as.numeric(dev200))), "\\n\\n")',
    '',
    '# --- Benchmarks ---',
    'b_cid_50 <- bench::mark(',
    '  ref = TreeDistRef::ClusteringInfoDistance(trees50),',
    '  dev = TreeDistDev::ClusteringInfoDistance(trees50),',
    '  min_iterations = 5, check = FALSE',
    ')',
    'cat("CID, 100 trees x 50 tips (4 950 pairs)\\n")',
    'print(b_cid_50[, c("expression", "min", "median", "mem_alloc", "n_itr")])',
    'cat("\\n")',
    '',
    'b_cid_200 <- bench::mark(',
    '  ref = TreeDistRef::ClusteringInfoDistance(trees200),',
    '  dev = TreeDistDev::ClusteringInfoDistance(trees200),',
    '  min_iterations = 5, check = FALSE',
    ')',
    'cat("CID, 40 trees x 200 tips (780 pairs)\\n")',
    'print(b_cid_200[, c("expression", "min", "median", "mem_alloc", "n_itr")])',
    'cat("\\n")',
    '',
    'b_msd_50 <- bench::mark(',
    '  ref = TreeDistRef::MatchingSplitDistance(trees50),',
    '  dev = TreeDistDev::MatchingSplitDistance(trees50),',
    '  min_iterations = 5, check = FALSE',
    ')',
    'cat("MSD, 100 trees x 50 tips (4 950 pairs)\\n")',
    'print(b_msd_50[, c("expression", "min", "median", "mem_alloc", "n_itr")])',
    'cat("\\n")',
    '',
    'b_msd_200 <- bench::mark(',
    '  ref = TreeDistRef::MatchingSplitDistance(trees200),',
    '  dev = TreeDistDev::MatchingSplitDistance(trees200),',
    '  min_iterations = 5, check = FALSE',
    ')',
    'cat("MSD, 40 trees x 200 tips (780 pairs)\\n")',
    'print(b_msd_200[, c("expression", "min", "median", "mem_alloc", "n_itr")])',
    'cat("\\n")',
    '',
    '# --- LAPJV standalone ---',
    'UnifMat <- function(n) matrix(runif(n * n), n, n)',
    'test400 <- UnifMat(400)',
    'test2000 <- UnifMat(1999)',
    '',
    'b_lap_400 <- bench::mark(',
    '  ref = TreeDistRef::LAPJV(test400),',
    '  dev = TreeDistDev::LAPJV(test400),',
    '  min_iterations = 5, check = FALSE',
    ')',
    'cat("LAPJV 400x400\\n")',
    'print(b_lap_400[, c("expression", "min", "median", "mem_alloc", "n_itr")])',
    'cat("\\n")',
    '',
    'b_lap_2000 <- bench::mark(',
    '  ref = TreeDistRef::LAPJV(test2000),',
    '  dev = TreeDistDev::LAPJV(test2000),',
    '  min_iterations = 5, check = FALSE',
    ')',
    'cat("LAPJV 1999x1999\\n")',
    'print(b_lap_2000[, c("expression", "min", "median", "mem_alloc", "n_itr")])',
    'cat("\\n")',
    '',
    sprintf('saveRDS(list(cid_50 = b_cid_50, cid_200 = b_cid_200, msd_50 = b_msd_50, msd_200 = b_msd_200, lap_400 = b_lap_400, lap_2000 = b_lap_2000), "%s")',
            gsub("\\\\", "/", out_rds)),
    sprintf('cat("Results saved to %s\\n")', gsub("\\\\", "/", out_rds))
  ), script)

  message("── Running A/B comparison in subprocess …")
  ret <- system2("Rscript", script)
  if (ret != 0L) stop("Benchmark subprocess failed (exit code ", ret, ")")

  invisible(readRDS(out_rds))
}
