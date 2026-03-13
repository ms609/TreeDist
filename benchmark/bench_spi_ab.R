# bench_spi_ab.R
# One-off A/B benchmark for the spi_overlap single-pass optimization.
# Uses the same ref lib as compare-ab.R; installs current dev to a temp lib.

source("benchmark/compare-ab.R")

# Install dev (current working tree with spi_overlap changes) to temp lib
stale <- Sys.glob(file.path(.REPO, "src", c("*.o", "*.dll", "*.so")))
if (length(stale)) file.remove(stale)

dev_lib <- file.path(tempdir(), "TreeDist_spi_dev")
if (!dir.exists(dev_lib)) dir.create(dev_lib, recursive = TRUE)

message("Installing dev TreeDist …")
install.packages(.REPO, lib = dev_lib, repos = NULL, type = "source",
                 INSTALL_opts = "--no-multiarch", quiet = FALSE)

script <- tempfile(fileext = ".R")
out_rds <- file.path(.RESDIR, "ab_spi.Rds")

writeLines(c(
  sprintf('.libPaths(c("%s", "%s", .libPaths()))',
          gsub("\\\\", "/", dev_lib),
          gsub("\\\\", "/", .REF_LIB)),
  'suppressPackageStartupMessages({',
  '  library(TreeTools)',
  '  library(TreeDistRef)',
  '  library(TreeDist)',
  '})',
  'cat("TreeDistRef:", find.package("TreeDistRef"), "\\n")',
  'cat("TreeDist   :", find.package("TreeDist"),    "\\n\\n")',
  '',
  'set.seed(4812)',
  'trees50  <- as.phylo(0:99, tipLabels = paste0("t", seq_len(50)))',
  'trees200 <- as.phylo(0:39, tipLabels = paste0("t", seq_len(200)))',
  '',
  '# Correctness checks',
  'chk <- function(ref, dev, label) {',
  '  d <- max(abs(as.numeric(ref) - as.numeric(dev)))',
  '  cat(sprintf("Max |ref-dev| %s: %.3e\\n", label, d))',
  '}',
  'chk(TreeDistRef::PhylogeneticInfoDistance(trees50),',
  '    TreeDist::PhylogeneticInfoDistance(trees50),  "PID  50-tip")',
  'chk(TreeDistRef::PhylogeneticInfoDistance(trees200),',
  '    TreeDist::PhylogeneticInfoDistance(trees200), "PID 200-tip")',
  'chk(TreeDistRef::ClusteringInfoDistance(trees50),',
  '    TreeDist::ClusteringInfoDistance(trees50),    "CID  50-tip (canary)")',
  'cat("\\n")',
  '',
  '# Benchmarks',
  'b_pid_50 <- bench::mark(',
  '  ref = TreeDistRef::PhylogeneticInfoDistance(trees50),',
  '  dev = TreeDist::PhylogeneticInfoDistance(trees50),',
  '  min_iterations = 5, check = FALSE',
  ')',
  'cat("PID, 100 trees x 50 tips (4 950 pairs)\\n")',
  'print(b_pid_50[, c("expression", "min", "median", "mem_alloc", "n_itr")])',
  'cat("\\n")',
  '',
  'b_pid_200 <- bench::mark(',
  '  ref = TreeDistRef::PhylogeneticInfoDistance(trees200),',
  '  dev = TreeDist::PhylogeneticInfoDistance(trees200),',
  '  min_iterations = 5, check = FALSE',
  ')',
  'cat("PID, 40 trees x 200 tips (780 pairs)\\n")',
  'print(b_pid_200[, c("expression", "min", "median", "mem_alloc", "n_itr")])',
  'cat("\\n")',
  '',
  '# CID as canary (should be unchanged)',
  'b_cid_50 <- bench::mark(',
  '  ref = TreeDistRef::ClusteringInfoDistance(trees50),',
  '  dev = TreeDist::ClusteringInfoDistance(trees50),',
  '  min_iterations = 5, check = FALSE',
  ')',
  'cat("CID, 100 trees x 50 tips (canary)\\n")',
  'print(b_cid_50[, c("expression", "min", "median", "mem_alloc", "n_itr")])',
  'cat("\\n")',
  '',
  sprintf('saveRDS(list(pid_50=b_pid_50, pid_200=b_pid_200, cid_50=b_cid_50), "%s")',
          gsub("\\\\", "/", out_rds)),
  sprintf('cat("Saved to %s\\n")', gsub("\\\\", "/", out_rds))
), script)

message("── Running A/B comparison in subprocess …")
ret <- system2("Rscript", script)
if (ret != 0L) stop("Subprocess failed (exit code ", ret, ")")

res <- readRDS(out_rds)
invisible(res)
