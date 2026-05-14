## Build shapeExpectedDistances lookup table
##
## For each n in 4:10, computes E[CID | shape1, shape2] and
## E[PID | shape1, shape2] over all (shape1 <= shape2) pairs by
## label-shuffle Monte Carlo.  Results are stored as a list keyed by n
## in inst/extdata/shapeExpectedDistances.rds.
##
## Shape convention (confirmed by round-trip): shapes are 0-based integers
## in 0:(NRootedShapes(n) - 1).  RootedTreeWithShape() requires Preorder()
## before passing to CID/PID functions.
##
## Usage:
##   Rscript data-raw/build_shape_lookup.R
##
## Wall-time budget: ~10 min per n at nSim = 1000.  Script prints cumulative
## file size after each n and stops if it would exceed 250 KB.

devtools::load_all("C:/Users/pjjg18/GitHub/worktrees/acid-aphid-TreeTools",
                   quiet = TRUE)
devtools::load_all(".", quiet = TRUE)

n_range  <- 4:10
nSim_default <- 1000L
size_cap_bytes <- 250000L          # 250 KB hard stop
out_dir  <- file.path("inst", "extdata")
out_file <- file.path(out_dir, "shapeExpectedDistances.rds")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------------
# Verify shape key convention
# -------------------------------------------------------------------------
cat("--- Shape key convention check ---\n")
for (n_chk in c(4L, 8L)) {
  W <- as.integer(TreeTools::NRootedShapes(n_chk))
  shapes_found <- vapply(seq_len(W) - 1L, function(s) {
    tr <- TreeTools::RootedTreeWithShape(s, n_chk)
    as.integer(TreeTools::RootedTreeShape(tr))
  }, integer(1L))
  cat(sprintf("  n=%d: NRootedShapes=%d, round-trip range=[%d,%d], ",
              n_chk, W, min(shapes_found), max(shapes_found)))
  if (isTRUE(all.equal(shapes_found, seq_len(W) - 1L))) {
    cat("0-based confirmed.\n")
  } else {
    cat("UNEXPECTED round-trip!\n")
    print(shapes_found)
  }
}
cat("\n")

# -------------------------------------------------------------------------
# Helper: MC estimate for one (s1, s2) pair
# -------------------------------------------------------------------------
.ShapePairMC <- function(s1, s2, n, nSim) {
  lab_pool <- as.character(seq_len(n))
  cid_v <- numeric(nSim)
  pid_v <- numeric(nSim)
  for (i in seq_len(nSim)) {
    t1 <- TreeTools::Preorder(
            TreeTools::RootedTreeWithShape(s1, n, sample(lab_pool)))
    t2 <- TreeTools::Preorder(
            TreeTools::RootedTreeWithShape(s2, n, sample(lab_pool)))
    cid_v[i] <- ClusteringInfoDistance(t1, t2, normalize = TRUE)
    pid_v[i] <- PhylogeneticInfoDistance(t1, t2, normalize = TRUE)
  }
  data.frame(
    shape1   = s1,
    shape2   = s2,
    mean_cid = mean(cid_v),
    sd_cid   = sd(cid_v),
    mean_pid = mean(pid_v),
    sd_pid   = sd(pid_v),
    n_pairs  = nSim
  )
}

# -------------------------------------------------------------------------
# Main loop
# -------------------------------------------------------------------------
result_list <- vector("list", length(n_range))
names(result_list) <- as.character(n_range)

t_total_start <- proc.time()[["elapsed"]]

for (n in n_range) {
  W <- as.integer(TreeTools::NRootedShapes(n))
  pairs_total <- W * (W + 1L) / 2L
  cat(sprintf("n=%d: W=%d shapes, %d pairs (nSim=%d per pair)\n",
              n, W, pairs_total, nSim_default))

  t_n_start <- proc.time()[["elapsed"]]

  rows <- vector("list", pairs_total)
  k <- 0L
  for (s1 in seq_len(W) - 1L) {
    for (s2 in s1:(W - 1L)) {
      k <- k + 1L
      rows[[k]] <- .ShapePairMC(s1, s2, n, nSim_default)
    }
    if (s1 %% 10L == 0L && s1 > 0L) {
      elapsed <- proc.time()[["elapsed"]] - t_n_start
      cat(sprintf("  s1=%d/%d, %d pairs done, %.1f s elapsed\n",
                  s1, W - 1L, k, elapsed))
    }
  }

  df <- do.call(rbind, rows)
  result_list[[as.character(n)]] <- df

  elapsed_n <- proc.time()[["elapsed"]] - t_n_start
  cat(sprintf("  n=%d done: %d rows in %.1f s\n", n, nrow(df), elapsed_n))

  # Write incrementally and check size
  saveRDS(result_list[!vapply(result_list, is.null, logical(1L))],
          out_file, compress = TRUE)
  sz <- file.size(out_file)
  cat(sprintf("  Cumulative file size: %.1f KB\n\n", sz / 1024))

  if (sz > size_cap_bytes) {
    cat(sprintf("!! File size %.1f KB exceeds cap of %.1f KB after n=%d.\n",
                sz / 1024, size_cap_bytes / 1024, n))
    cat("   Stopping here; shipped table covers n=",
        paste(n_range[n_range <= n], collapse = ":"), "\n")
    break
  }
}

# -------------------------------------------------------------------------
# Final summary
# -------------------------------------------------------------------------
final_list <- result_list[!vapply(result_list, is.null, logical(1L))]
saveRDS(final_list, out_file, compress = TRUE)

total_elapsed <- proc.time()[["elapsed"]] - t_total_start
cat(sprintf("=== Done: n covered = %s, file = %.1f KB, wall time = %.1f s ===\n",
            paste(names(final_list), collapse = ":"),
            file.size(out_file) / 1024,
            total_elapsed))
