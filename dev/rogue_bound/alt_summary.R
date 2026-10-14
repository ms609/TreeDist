dir <- "/tmp/claude-0/-home-claude/37427f5a-6383-50c4-a6e2-4664447d3c94/scratchpad/bound/csv2"
out <- NULL
for (f in c("rand30", "rand30d2", "rand80", "rand150", "near80", "near150")) {
  x <- read.csv(file.path(dir, paste0(f, ".csv")))
  lapCalls <- x$ub_row >= x$min_score
  w <- x[lapCalls & x$lap_exact < x$min_score, ]
  Rej <- function(b) { v <- w[[b]]; v[is.na(v)] <- Inf; 100 * mean(v < w$min_score) }
  out <- rbind(out, data.frame(
    wl = f, n_wasted = nrow(w),
    cr = Rej("ub_cr"), rc = Rej("ub_rc"), best_cr_rc = Rej("ub_alt1"),
    alt3 = Rej("ub_alt3"), alt10 = Rej("ub_alt10"),
    pv = Rej("ub_pv_minv"), pv_max = Rej("ub_pv_maxv"), pv_mean = Rej("ub_pv_mean"),
    pu = Rej("ub_pu"), pv1 = Rej("ub_pv_alt1"), pv3 = Rej("ub_pv_alt3"),
    pvarr = Rej("ub_pv_warmarr"), prev = Rej("ub_prev"),
    arr = Rej("ub_arr"), rt = Rej("ub_rt")))
}
print(format(out, digits = 3), row.names = FALSE)
