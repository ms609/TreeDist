# Usage: Rscript verify.R <lib> <out.rds> [capture.rds]
# Runs SPRDist(method = "Rogue") on the verify pairs at depths 1 and 2.
# With a third argument, also records every spr_rogue_scores() input
# (from this lib's run) to that file.
args <- commandArgs(TRUE)
.libPaths(c("/root/Rlib", .libPaths()))
library(TreeDist, lib.loc = args[1])
suppressMessages(library(TreeTools))
pairs <- readRDS("/tmp/claude-0/-home-claude/37427f5a-6383-50c4-a6e2-4664447d3c94/scratchpad/bound/verify.rds")
capture <- length(args) > 2
captured <- list()
if (capture) {
  ns <- asNamespace("TreeDist")
  orig <- get("spr_rogue_scores", ns)
  env <- environment()
  wrapped <- function(tree1, tree2, depth) {
    assign("captured", c(get("captured", env), list(list(tree1, tree2, depth))), env)
    orig(tree1, tree2, depth)
  }
  unlockBinding("spr_rogue_scores", ns)
  assign("spr_rogue_scores", wrapped, ns)
}
nTip <- vapply(pairs, function(p) NTip(p[[1]]), 1)
options(sprDepth = 1)
d1 <- vapply(pairs, function(p) SPRDist(p[[1]], p[[2]], method = "Rogue"), 1)
options(sprDepth = 2)
small <- which(nTip <= 60)
d2 <- vapply(pairs[small], function(p) SPRDist(p[[1]], p[[2]], method = "Rogue"), 1)
saveRDS(list(d1 = d1, d2 = d2), args[2])
if (capture) saveRDS(captured, args[3])
cat(basename(args[1]), "pairs:", length(d1), "depth-2 pairs:", length(d2),
    "sum d1", sum(d1), "sum d2", sum(d2), "captured", length(captured), "\n")
