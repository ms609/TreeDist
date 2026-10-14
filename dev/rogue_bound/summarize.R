# Summarize analyzer CSVs: which bounds would reject calls that currently
# need a LAP only to be rejected.
args <- commandArgs(TRUE)
dir <- "/tmp/claude-0/-home-claude/37427f5a-6383-50c4-a6e2-4664447d3c94/scratchpad/bound/csv"
files <- if (length(args)) args else sub("\\.csv$", "", list.files(dir, "\\.csv$"))
bounds <- c("ub_entropy", "ub_row", "ub_cr", "ub_rc", "ub_rt", "ub_arr",
            "ub_alt1", "ub_alt3", "ub_alt10", "ub_pv_minv", "ub_pv_maxv",
            "ub_pv_mean", "ub_pu", "ub_pv_alt1", "ub_pv_alt3", "ub_pv_warmarr",
            "ub_prev", "ub_size")
for (f in files) {
  x <- read.csv(file.path(dir, paste0(f, ".csv")))
  x$ub_size[is.na(x$ub_size)] <- Inf
  for (b in bounds) x[[b]][is.na(x[[b]])] <- Inf
  # Sanity: every bound must be >= the exact score
  for (b in bounds) {
    bad <- sum(x[[b]] < x$lap_exact - 1e-9)
    if (bad) cat("!! INVALID", b, bad, "\n")
  }
  lapCalls <- x$ub_row >= x$min_score  # not caught by the current row bound
  wasted <- lapCalls & x$lap_exact < x$min_score
  cat(sprintf("\n=== %s: %d calls; reach LAP %d; LAP only rejects %d (%.1f%%)\n",
              f, nrow(x), sum(lapCalls), sum(wasted), 100 * mean(wasted[lapCalls])))
  w <- x[wasted, ]
  rel <- (w$min_score - w$lap_exact) / w$lap_exact
  cat(sprintf("  margin (min_score - MCI) for wasted LAPs: median %.4f bits (%.3f%% of MCI); q10 %.4f, q90 %.4f\n",
              median(w$min_score - w$lap_exact), 100 * median(rel),
              quantile(w$min_score - w$lap_exact, 0.1), quantile(w$min_score - w$lap_exact, 0.9)))
  res <- t(vapply(bounds, function(b) {
    c(rejects = mean(w[[b]] < w$min_score),
      medExcess = median(w[[b]] - w$lap_exact),
      medExcessPct = 100 * median((w[[b]] - w$lap_exact) / w$lap_exact))
  }, double(3)))
  print(round(res, 4))
  if (nrow(w)) {
    cat(sprintf("  augmentation phase: median free rows after ARR %.0f; bound rejects within aug %.1f%% of wasted (median frac of augs needed %.2f)\n",
                median(w$aug_total), 100 * mean(w$aug_needed >= 0 & w$aug_needed < w$aug_total),
                median(ifelse(w$aug_needed >= 0, w$aug_needed / pmax(1, w$aug_total), 1))))
    cat(sprintf("  frac of cols mapped from parent: median %.2f; multi-preimage cols median %.0f\n",
                median(x$frac_mapped, na.rm = TRUE), median(x$n_multi, na.rm = TRUE)))
    lx <- x[lapCalls, ]
    cat(sprintf("  LAP work (all LAP calls): cold ticks %.0f, warm ticks %.0f (ratio %.2f); warm optimal %.4f\n",
                mean(lx$cold_ticks), mean(lx$warm_ticks, na.rm = TRUE),
                sum(lx$warm_ticks, na.rm = TRUE) / sum(lx$cold_ticks[!is.na(lx$warm_ticks)]),
                mean(lx$warm_ok[!is.na(lx$warm_ticks) & lx$warm_ticks > 0])))
    cat(sprintf("  cold: free rows after init %.1f, after ARR %.1f, aug ops %.0f; warm: %.1f, %.1f, %.0f\n",
                mean(lx$cold_free_init), mean(lx$cold_free_arr), mean(lx$cold_aug),
                mean(lx$warm_free_init), mean(lx$warm_free_arr), mean(lx$warm_aug)))
    cat(sprintf("  mean d %.1f\n", mean(lx$d)))
  }
}
