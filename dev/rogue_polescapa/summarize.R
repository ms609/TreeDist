# Tables for dev/rogue_polescapa_findings.md
# Rscript dev/rogue_polescapa/summarize.R
sel <- c("mci", "both-k0", "both-k0-last", "12-k0", "21-k0", "both-k3",
         "both-k5", "both-k10", "both-tied", "12-tied", "mci-then-fitch",
         "hyb20-k0", "hyb20-k5")

# Results, with the hybrid columns from the `hyb` runs merged in
Load <- function(which) {
  f <- sprintf("dev/rogue_polescapa/results-%s.rds", which)
  if (!file.exists(f)) return(NULL)
  out <- readRDS(f)
  h <- sprintf("dev/rogue_polescapa/results-%s-hyb.rds", which)
  if (file.exists(h)) {
    hyb <- readRDS(h)$res
    stopifnot(identical(hyb$mci, out$res$mci))
    out$res[, c("hyb20-k0", "hyb20-k5")] <- hyb[, c("hyb20-k0", "hyb20-k5")]
  }
  out
}

Dist <- function(x) {
  x <- x[!is.na(x)]
  sprintf("%d/%d/%d/%d (%.2f)", sum(x == 0), sum(x == 1), sum(x == 2),
          sum(x >= 3), mean(x))
}

Table <- function(res, refCol, title) {
  cat("\n###", title, "\n\n")
  cat("Cells: pairs with over-estimate 0 / +1 / +2 / +3 or more (mean).",
      "W/L: pairs where the selector beats / loses to `mci`.\n\n")
  sets <- unique(res$set)
  cat("| Selector |", paste(sets, collapse = " | "), "| W/L vs mci |\n")
  cat("|", paste(rep("---", length(sets) + 2), collapse = " | "), "|\n")
  for (s in intersect(sel, names(res))) {
    cells <- vapply(sets, function(st) {
      r <- res[res$set == st, ]
      Dist(r[[s]] - r[[refCol]])
    }, "")
    wl <- sprintf("%d/%d", sum(res[[s]] < res$mci), sum(res[[s]] > res$mci))
    cat("|", s, "|", paste(cells, collapse = " | "), "|", wl, "|\n")
  }
  cat("\nPairs per set:", paste(sets, table(res$set)[sets], collapse = ", "),
      "\n")
}

TieTable <- function(log, title) {
  cat("\n###", title, "\n\n")
  lg <- log[log$selector == "both-k0", ]
  cat("| Set | Iterations | Median n | Median tie set | Tie set > 1 | Tie set >= 5 |\n")
  cat("| --- | --- | --- | --- | --- | --- |\n")
  for (st in unique(lg$set)) {
    l <- lg[lg$set == st, ]
    cat(sprintf("| %s | %d | %d | %d | %.0f%% | %.0f%% |\n", st, nrow(l),
                as.integer(median(l$n)), as.integer(median(l$nTied)),
                100 * mean(l$nTied > 1), 100 * mean(l$nTied >= 5)))
  }
}

ref <- Load("reference")
r <- ref$res
cat("Reference pairs with exact distance:", sum(!is.na(r$exact)), "of",
    nrow(r), "\n")
Table(r[!is.na(r$exact), ], "exact", "Over-estimate vs exact uSPR (TBRDist)")
TieTable(ref$log, "Fitch tie sets (both directions), reference pairs")

for (which in c("large", "extra")) {
  lg <- Load(which)
  if (is.null(lg)) next
  l <- lg$res
  l$best <- apply(l[, intersect(sel, names(l))], 1, min)
  Table(l, "best", paste0("Excess over best-known (min over selectors), ",
                          which, " set"))
  TieTable(lg$log, paste("Fitch tie sets (both directions),", which, "set"))
}
