# SPR estimates under each selector.
# Usage (from package root):
#   Rscript dev/rogue_polescapa/evaluate.R reference   # pairs with exact uSPR
#   Rscript dev/rogue_polescapa/evaluate.R large       # dev/rogue_bound workloads
Sys.setenv(MAKEFLAGS = "-j2")
suppressMessages(devtools::load_all(quiet = TRUE))
source("dev/rogue_polescapa/selectors.R")
#   Rscript dev/rogue_polescapa/evaluate.R extra       # more large random pairs
# A second argument `hyb` runs only `mci` and the hybrids, saved separately.
args <- commandArgs(TRUE)
which <- args[[1]]
hybOnly <- identical(args[2], "hyb")

selectors <- list(
  mci = ScorerMCI,
  `both-k0` = MakeFitchScorer("both", 0),
  `both-k0-last` = MakeFitchScorer("both", 0, last = TRUE),
  `12-k0` = MakeFitchScorer("12", 0),
  `21-k0` = MakeFitchScorer("21", 0),
  `both-k3` = MakeFitchScorer("both", 3),
  `both-k5` = MakeFitchScorer("both", 5),
  `both-k10` = MakeFitchScorer("both", 10),
  `both-tied` = MakeFitchScorer("both", "tied"),
  `12-tied` = MakeFitchScorer("12", "tied"),
  `mci-then-fitch` = ScorerMCIThenFitch,
  `hyb20-k0` = MakeHybridScorer(20, 0),
  `hyb20-k5` = MakeHybridScorer(20, 5)
)
if (hybOnly) {
  selectors <- selectors[c("mci", "hyb20-k0", "hyb20-k5")]
} else if (which == "extra") {
  selectors <- selectors[c("mci", "both-k0", "both-k5", "both-tied",
                           "hyb20-k0", "hyb20-k5")]
}

if (which == "reference") {
  ref <- readRDS("dev/rogue_polescapa/reference.rds")
  pairs <- lapply(ref, function(r) r[c("trees", "set", "nTip", "moves", "exact")])
} else if (which == "extra") {
  # Further large random pairs (seeded), to firm up the 80/150-leaf results
  source("dev/rogue_polescapa/workloads.R")
  set.seed(20261024)
  Make <- function(nPair, nTip, set) lapply(seq_len(nPair), function(i) {
    list(trees = RandomPair(nTip), set = set, nTip = nTip, moves = NA,
         exact = NA)
  })
  pairs <- c(Make(30, 50, "rand50"), Make(20, 80, "rand80"),
             Make(10, 150, "rand150"))
} else {
  wl <- readRDS("dev/rogue_bound/workloads.rds")
  pairs <- unlist(lapply(names(wl), function(set) {
    lapply(wl[[set]], function(tr) {
      list(trees = tr, set = set, nTip = NTip(tr[[1]]), moves = NA,
           exact = NA)
    })
  }), recursive = FALSE)
}

rows <- list()
logs <- list()
for (i in seq_along(pairs)) {
  pr <- pairs[[i]]
  est <- vapply(names(selectors), function(s) {
    ResetLog()
    t0 <- proc.time()[["elapsed"]]
    d <- RogueWith(pr$trees[[1]], pr$trees[[2]], selectors[[s]])
    el <- proc.time()[["elapsed"]] - t0
    lg <- .log$rows
    if (length(lg)) {
      logs[[length(logs) + 1L]] <<- data.frame(
        pair = i, set = pr$set, selector = s,
        n = vapply(lg, `[[`, 0, "n"),
        nTied = vapply(lg, function(r) if (is.null(r$nTied)) NA else r$nTied, 0)
      )
    }
    d
  }, double(1))
  rows[[i]] <- data.frame(pair = i, set = pr$set, nTip = pr$nTip,
                          moves = pr$moves, exact = pr$exact, t(est),
                          check.names = FALSE)
  if (i %% 25 == 0) cat(i, "/", length(pairs), "\n")
}
res <- do.call(rbind, rows)
saveRDS(list(res = res, log = do.call(rbind, logs)),
        sprintf("dev/rogue_polescapa/results-%s%s.rds", which,
                if (hybOnly) "-hyb" else ""))
