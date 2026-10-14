# Usage: Rscript scores.R <lib> <capture.rds> <out.rds>
# Calls spr_rogue_scores() on each captured input; saves the outputs.
args <- commandArgs(TRUE)
.libPaths(c("/root/Rlib", .libPaths()))
library(TreeDist, lib.loc = args[1])
inputs <- readRDS(args[2])
out <- lapply(inputs, function(x) TreeDist:::spr_rogue_scores(x[[1]], x[[2]], x[[3]]))
saveRDS(out, args[3])
cat(basename(args[1]), length(out), "calls\n")
