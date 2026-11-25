# Checks (1) the incremental Fitch against TreeSearch::TreeLength() and
# (2) that RogueWith(ScorerMCI) reproduces SPRDist(method = "Rogue").
# Run from the package root: Rscript dev/rogue_polescapa/validate.R
Sys.setenv(MAKEFLAGS = "-j2")
suppressMessages(devtools::load_all(quiet = TRUE))
source("dev/rogue_polescapa/selectors.R")
source("dev/rogue_polescapa/workloads.R")
suppressMessages(library("TreeSearch"))

# (1) Fitch: length with each leaf coded `?` in every character
FitchByTreeLength <- function(tree, chars) {
  mat <- t(ifelse(chars, "1", "0"))
  rownames(mat) <- TipLabels(tree)
  tree <- RootTree(tree, 1)
  dat <- MatrixToPhyDat(mat)
  full <- TreeLength(tree, dat)
  amb <- vapply(seq_len(NTip(tree)), function(x) {
    m <- mat
    m[x, ] <- "?"
    TreeLength(tree, MatrixToPhyDat(m))
  }, double(1))
  list(length = full, ambiguous = amb)
}

set.seed(1)
fitchOK <- vapply(1:40, function(i) {
  n <- sample(8:40, 1)
  pr <- if (i %% 2) RandomPair(n) else NearPair(n, sample(1:4, 1))
  chars <- SplitChars(pr[[1]], pr[[2]])
  a <- FitchAmbiguous(pr[[1]], chars)
  b <- FitchByTreeLength(pr[[1]], chars)
  isTRUE(all.equal(a$length, b$length)) &&
    isTRUE(all.equal(as.numeric(a$ambiguous), as.numeric(b$ambiguous)))
}, logical(1))
cat("Incremental Fitch == TreeLength:", sum(fitchOK), "/", length(fitchOK), "\n")

# (2) Harness reproduces SPRDist()
set.seed(2)
pairs <- c(lapply(1:30, function(i) RandomPair(sample(8:60, 1))),
           lapply(1:30, function(i) NearPair(sample(10:80, 1), sample(1:8, 1))))
same <- vapply(pairs, function(pr) {
  SPRDist(pr[[1]], pr[[2]], method = "Rogue")[[1]] ==
    RogueWith(pr[[1]], pr[[2]], ScorerMCI)
}, logical(1))
cat("Harness == SPRDist(Rogue):", sum(same), "/", length(same), "\n")
