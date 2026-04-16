#!/usr/bin/env Rscript
# Quick test runner for TreeDist

library("devtools", quietly = TRUE)
library("testthat", quietly = TRUE)

cat("Building and loading TreeDist...\n")
load_all()

cat("\nRunning test suite (quick mode, no vignettes)...\n")
test()
