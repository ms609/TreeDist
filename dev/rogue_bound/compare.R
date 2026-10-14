# Usage: Rscript compare.R <base_scores.rds> <new_scores.rds>
# Checks, for each spr_rogue_scores() call, that both versions agree on the
# minimum score and on the set of candidates attaining it, and that every
# candidate scored exactly by either version (score within tolerance of the
# minimum) has the same score in both.
args <- commandArgs(TRUE)
a <- readRDS(args[1])
b <- readRDS(args[2])
stopifnot(length(a) == length(b))
Check <- function(x, y, label) {
  if (length(x) != length(y)) return(sprintf("%s: length", label))
  if (!length(x)) return(NULL)
  fin <- is.finite(x) & is.finite(y)
  if (!identical(is.finite(x), is.finite(y))) return(sprintf("%s: finiteness", label))
  if (!any(fin)) return(NULL)
  mx <- min(x[fin]); my <- min(y[fin])
  tol <- 1e-9 * max(1, abs(mx))
  if (abs(mx - my) > tol) return(sprintf("%s: min %.15g vs %.15g", label, mx, my))
  # Candidates tied with the minimum (exact equality, as SPRDist uses)
  if (!identical(which(x == mx), which(y == my))) {
    near <- abs(x - mx) <= tol | abs(y - my) <= tol
    return(sprintf("%s: argmin set differs (%s vs %s); near-ties: %d",
                   label, paste(which(x == mx), collapse = " "),
                   paste(which(y == my), collapse = " "), sum(near)))
  }
  # Every candidate close to the minimum must be exact in both
  close <- fin & (x <= mx + 0.5 | y <= my + 0.5)
  both <- close & x <= mx + 1e-6 & y <= my + 1e-6
  if (any(abs(x[both] - y[both]) > tol)) return(sprintf("%s: near-min scores differ", label))
  # A candidate rejected by one version must be rejected (> min) by the other
  if (any(x[fin] <= mx & y[fin] > my) || any(y[fin] <= my & x[fin] > mx)) {
    return(sprintf("%s: rejection mismatch", label))
  }
  NULL
}
problems <- character(0)
nCand <- 0
for (i in seq_along(a)) {
  p1 <- Check(a[[i]]$scores, b[[i]]$scores, sprintf("call %d singles", i))
  p2 <- Check(a[[i]]$pair_scores, b[[i]]$pair_scores, sprintf("call %d pairs", i))
  nCand <- nCand + length(a[[i]]$scores) + length(a[[i]]$pair_scores)
  problems <- c(problems, p1, p2)
}
cat("calls:", length(a), " candidates:", nCand, " problems:", length(problems), "\n")
if (length(problems)) print(head(problems, 20))
# How often exact scores differ at all (machine-precision level)
diffs <- unlist(lapply(seq_along(a), function(i) {
  x <- a[[i]]$scores; y <- b[[i]]$scores
  m <- min(x); keep <- is.finite(x) & x == m & y == min(y)
  abs(x[keep] - y[keep])
}))
cat("max |diff| among minimum scores:", max(c(0, diffs)), "\n")
