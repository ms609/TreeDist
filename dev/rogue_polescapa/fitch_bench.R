# Cost of the Fitch selector in C++ vs current per-iteration MCI scoring.
# Times one bit-parallel Fitch scoring (both directions, including the
# split extraction) against one spr_rogue_scores() call, at 30/80/150 leaves,
# and whole SPRDist(method = "Rogue") at HEAD.
# Run from the package root: Rscript dev/rogue_polescapa/fitch_bench.R
Sys.setenv(MAKEFLAGS = "-j2")
suppressMessages(devtools::load_all(quiet = TRUE))
source("dev/rogue_polescapa/selectors.R")
Rcpp::sourceCpp("dev/rogue_polescapa/fitch_bitpar.cpp")

FitchCpp <- function(tr1, tr2) {
  a <- Preorder(RootTree(tr1, 1))
  b <- Preorder(RootTree(tr2, 1))
  fitch_added(a[["edge"]], NTip(a), SplitChars(a, b)) +
    fitch_added(b[["edge"]], NTip(b), SplitChars(b, a))[
      match(TipLabels(a), TipLabels(b))]
}

# Check against the R prototype
wl <- readRDS("dev/rogue_bound/workloads.rds")
ok <- vapply(wl$rand30, function(pr) {
  a <- Preorder(RootTree(pr[[1]], 1))
  b <- Preorder(RootTree(pr[[2]], 1))
  all(FitchCpp(a, b) == FitchScores(a, b, "both"))
}, logical(1))
cat("C++ == R prototype:", sum(ok), "/", length(ok), "\n")

Med <- function(expr, reps) {
  e <- substitute(expr)
  env <- parent.frame()
  t <- replicate(5, system.time(
    for (i in seq_len(reps)) eval(e, env))[["elapsed"]])
  median(t) / reps
}
for (set in c("rand30", "rand80", "rand150")) {
  pr <- wl[[set]][[1]]
  red <- ReduceTrees(pr[[1]], pr[[2]])
  a <- red[[1]]; b <- red[[2]]
  aa <- Preorder(RootTree(a, 1)); bb <- Preorder(RootTree(b, 1))
  ca <- SplitChars(aa, bb); cb <- SplitChars(bb, aa)
  reps <- 2000
  kernel <- Med({fitch_added(aa[["edge"]], NTip(aa), ca)
                 fitch_added(bb[["edge"]], NTip(bb), cb)}, reps)
  withR <- Med(FitchCpp(a, b), 200)
  mci <- Med(spr_rogue_scores(a, b, 1L), 20)
  whole <- Med(vapply(wl[[set]], function(p) SPRDist(p[[1]], p[[2]],
                                                    method = "Rogue")[[1]], 0),
               1) / length(wl[[set]])
  cat(sprintf(
    "%s: n=%d  Fitch kernel %.1f us | Fitch incl. R glue %.0f us | spr_rogue_scores %.0f us | SPRDist/pair %.1f ms\n",
    set, NTip(a), kernel * 1e6, withR * 1e6, mci * 1e6, whole * 1e3))
}

# End to end: harness with the C++ Fitch selector vs SPRDist() at HEAD.
# The Fitch figure includes the R glue (as.Splits, RootTree) per iteration.
ScorerFitchCpp <- function(tr1, tr2) which.max(FitchCpp(tr1, tr2))
for (set in c("rand30", "rand80", "rand150", "mid80", "mid150")) {
  prs <- wl[[set]]
  tF <- system.time(f <- vapply(prs, function(p)
    RogueWith(p[[1]], p[[2]], ScorerFitchCpp), 0))[["elapsed"]]
  tM <- system.time(m <- vapply(prs, function(p)
    SPRDist(p[[1]], p[[2]], method = "Rogue")[[1]], 0))[["elapsed"]]
  cat(sprintf("%s: per pair  Fitch harness %.1f ms | SPRDist %.1f ms | mean est %.2f vs %.2f\n",
              set, 1e3 * tF / length(prs), 1e3 * tM / length(prs),
              mean(f), mean(m)))
}
