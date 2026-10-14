# Min time per (lib, workload, depth) across interleaved rounds; % vs base.
x <- read.table(text = grep("^lib_", readLines(
  "/tmp/claude-0/-home-claude/37427f5a-6383-50c4-a6e2-4664447d3c94/scratchpad/bound/bench_results.txt"),
  value = TRUE), stringsAsFactors = FALSE)
names(x) <- c("lib", "wl", "depth", "minLab", "min", "medLab", "med", "sum")
x$key <- paste(x$wl, x$depth)
mins <- tapply(x$min, list(x$key, x$lib), min)
runs <- tapply(x$min, list(x$key, x$lib), length)
out <- data.frame(base = mins[, "lib_base"], arr = mins[, "lib_arr"],
                  new2 = mins[, "lib_new2"],
                  arrPct = round(100 * (mins[, "lib_arr"] / mins[, "lib_base"] - 1)),
                  new2Pct = round(100 * (mins[, "lib_new2"] / mins[, "lib_base"] - 1)),
                  rounds = runs[, "lib_base"])
print(out)
sums <- tapply(x$sum, x$key, function(s) length(unique(s)))
cat("checksums agree across libs:", all(sums == 1), "\n")
