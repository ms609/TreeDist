# Usage: Rscript profile2.R <lib> <workload> [depth] [subset]
# For the instrumented prototype (td_new_instr).
args <- commandArgs(TRUE)
.libPaths(c("/root/Rlib", .libPaths()))
library(TreeDist, lib.loc = args[1])
suppressMessages(library(TreeTools))
wl <- readRDS("/tmp/claude-0/-home-claude/37427f5a-6383-50c4-a6e2-4664447d3c94/scratchpad/bound/workloads.rds")
pairs <- wl[[args[2]]]
if (length(args) > 3) pairs <- pairs[eval(parse(text = args[4]))]
depth <- if (length(args) > 2) as.integer(args[3]) else 1L
options(sprDepth = depth)
Run <- function() vapply(pairs, function(p) SPRDist(p[[1]], p[[2]], method = "Rogue"), 1)
invisible(Run())
invisible(TreeDist:::instr_get(TRUE))
el <- system.time(r <- Run())[["elapsed"]]
x <- TreeDist:::instr_get(TRUE)
timers <- c("total", "parentPrices", "childPrices", "keepReduce", "splits",
            "firstMatch", "mciFill", "mciLap", "boundedCID")
counts <- c("calls", "rowExit", "pricedAvail", "pricedExit", "lap", "lapBelow", "exactAll")
tk <- x[seq_along(timers)] / 2.1e9
names(tk) <- timers
ct <- x[length(timers) + seq_along(counts)]
names(ct) <- counts
cat(sprintf("== %s %s depth %d: elapsed %.3fs, sum %g\n", basename(args[1]), args[2], depth, el, sum(r)))
print(rbind(seconds = round(tk, 3), share = round(tk / el, 3)))
print(ct)
