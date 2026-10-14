# Progressive (row-by-row) priced bound: how much of the fill could be
# skipped if the parent-price bound were checked while rows are filled?
dir <- "/tmp/claude-0/-home-claude/37427f5a-6383-50c4-a6e2-4664447d3c94/scratchpad/bound/csv2"
for (f in c("rand30", "rand30d2", "rand80", "rand150", "near80", "near150")) {
  x <- read.csv(file.path(dir, paste0(f, ".csv")))
  lapCalls <- x$ub_row >= x$min_score
  wasted <- lapCalls & x$lap_exact < x$min_score
  priced <- lapCalls & !is.na(x$ub_pv_minv) & x$ub_pv_minv < x$min_score
  w <- x[priced, ]
  Frac <- function(p) ifelse(is.na(p) | p > 1, 1, p)
  # Fill cost ~ d^2; a reject after fraction f of rows saves (1 - f) of it.
  Saved <- function(col, rows) {
    sum((1 - Frac(x[[col]][rows])) * x$d[rows]^2) / sum(x$d[lapCalls]^2)
  }
  cat(sprintf("=== %s: LAP-stage calls %d, wasted %d, priced-bound rejects %d (%.1f%% of wasted)\n",
              f, sum(lapCalls), sum(wasted), sum(priced), 100 * sum(priced) / max(1, sum(wasted))))
  cat(sprintf("  rejected before any row is filled (size profile + parent prices): %.1f%% of wasted\n",
              100 * sum(lapCalls & !is.na(x$prog_nat) & x$prog_nat == 0) / max(1, sum(wasted))))
  for (k in c("prog_nat", "prog_h", "prog_oracle")) {
    pk <- w[[k]]
    cat(sprintf("  %-11s rows needed among priced rejects: median %.2f, q25 %.2f, q75 %.2f; fill saved (d^2-weighted, all LAP-stage calls) %.1f%%\n",
                k, median(Frac(pk)), quantile(Frac(pk), 0.25), quantile(Frac(pk), 0.75),
                100 * Saved(k, which(priced))))
  }
}
