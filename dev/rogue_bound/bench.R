# Usage: Rscript bench.R <lib> <workload> <reps> [depth]
# Prints min and median elapsed time over `reps` runs, and a checksum.
args <- commandArgs(TRUE)
.libPaths(c("/root/Rlib", .libPaths()))
library(TreeDist, lib.loc = args[1])
suppressMessages(library(TreeTools))
wl <- readRDS("/tmp/claude-0/-home-claude/37427f5a-6383-50c4-a6e2-4664447d3c94/scratchpad/bound/workloads.rds")
pairs <- wl[[args[2]]]
reps <- as.integer(args[3])
depth <- if (length(args) > 3) as.integer(args[4]) else 1L
options(sprDepth = depth)
Run <- function() vapply(pairs, function(p) SPRDist(p[[1]], p[[2]], method = "Rogue"), 1)
r <- Run()
el <- replicate(reps, system.time(Run())[["elapsed"]])
cat(sprintf("%-10s %-8s d=%d min %.3f med %.3f sum=%g\n", basename(args[1]),
            args[2], depth, min(el), median(el), sum(r)))
