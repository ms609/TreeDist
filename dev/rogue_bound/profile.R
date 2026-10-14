# Usage: Rscript profile.R <lib> <workload> [depth] [dumpfile]
args <- commandArgs(TRUE)
.libPaths(c("/root/Rlib", .libPaths()))
library(TreeDist, lib.loc = args[1])
suppressMessages(library(TreeTools))
wl <- readRDS("/tmp/claude-0/-home-claude/37427f5a-6383-50c4-a6e2-4664447d3c94/scratchpad/bound/workloads.rds")
pairs <- wl[[args[2]]]
if (length(args) > 4) pairs <- pairs[eval(parse(text = args[5]))]
depth <- if (length(args) > 2) as.integer(args[3]) else 1L
options(sprDepth = depth)
dumpFile <- if (length(args) > 3) args[4] else ""
Run <- function() vapply(pairs, function(p) SPRDist(p[[1]], p[[2]], method = "Rogue"), 1)
invisible(Run()) # warm up
invisible(TreeDist:::instr_get(TRUE))
if (nzchar(dumpFile)) TreeDist:::instr_dump(dumpFile)
el <- system.time(r <- Run())[["elapsed"]]
TreeDist:::instr_dump("")
x <- TreeDist:::instr_get(TRUE)
timers <- c("total", "keepReduce", "splits", "firstMatch", "boundedCID",
            "mciSetup", "mciFill", "mciCopy", "mciLap", "mciTail", "halfKeepReduce", "parent")
counts <- c("calls", "entropyExit", "rowExit", "exactAll", "lap", "lapBelow",
            "whole", "half1", "half2", "candidates", "extraExit")
tk <- x[seq_along(timers)] / 2.1e9
names(tk) <- timers
ct <- x[length(timers) + seq_along(counts)]
names(ct) <- counts
cat(sprintf("== %s %s depth %d: elapsed %.3fs, sum %g\n", basename(args[1]), args[2], depth, el, sum(r)))
cat("Seconds (C++ timers):\n"); print(round(tk, 3))
cat("Share of elapsed:\n"); print(round(tk / el, 3))
cat("Counts:\n"); print(ct)
