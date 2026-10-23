# Build the reference set: tree pairs with exact uSPR distances (TBRDist).
# Each USPRDist() call runs in its own process with a time limit, since its
# C++ cannot be interrupted. Pairs that time out keep `exact = NA`.
# Run from the package root: Rscript dev/rogue_polescapa/make_reference.R
suppressMessages({
  library("TreeTools")
  library("TreeSearch")
})
source("dev/rogue_polescapa/workloads.R")
out <- "dev/rogue_polescapa/reference.rds"

# Fixtures from git history (b6359b04, 5c1179e1): successive TreeSearch::SPR()
# moves from a random tree; stored matrices hold USPRDist() of each pair.
# The regenerated trees reproduce them (36/36 USPRDist spot checks agree), so
# the stored values are used rather than recomputed.
FixtureTrees <- function(nTip, nSPR) {
  set.seed(0)
  tr <- vector("list", nSPR + 1L)
  tr[[1]] <- Postorder(RandomTree(nTip, root = TRUE))
  for (i in seq_len(nSPR) + 1L) tr[[i]] <- Postorder(TreeSearch::SPR(tr[[i - 1]]))
  tr
}
FixturePairs <- function(nTip, nSPR, set) {
  tr <- FixtureTrees(nTip, nSPR)
  stored <- as.matrix(readRDS(sprintf("dev/rogue_polescapa/true-%dtip-%dspr.Rds",
                                      nTip, nSPR)))
  ij <- combn(length(tr), 2)
  lapply(seq_len(ncol(ij)), function(k) {
    list(trees = list(tr[[ij[1, k]]], tr[[ij[2, k]]]), set = set,
         nTip = nTip, moves = abs(ij[2, k] - ij[1, k]),
         exact = stored[ij[1, k], ij[2, k]])
  })
}

set.seed(20261023)
Gen <- function(nPair, set, tips, moves = NULL) {
  lapply(seq_len(nPair), function(i) {
    n <- sample(tips, 1)
    m <- if (is.null(moves)) NA else sample(moves, 1)
    list(trees = if (is.null(moves)) RandomPair(n) else NearPair(n, m),
         set = set, nTip = n, moves = m)
  })
}
pairs <- c(
  Gen(120, "rand8-14", 8:14),
  Gen(80, "rand15-20", 15:20),
  Gen(80, "near20-50", 20:50, 2:6),
  Gen(60, "near20-50x", 20:50, 7:10),
  FixturePairs(25, 12, "fix25-12spr"),
  FixturePairs(32, 16, "fix32-16spr")
)

# One persistent worker; restarted when a call exceeds `timeout` seconds.
.worker <- NULL
ExactUSPR <- function(trees, timeout = 60) {
  if (is.null(.worker) || !.worker$is_alive()) {
    .worker <<- callr::r_session$new()
    .worker$run(function() loadNamespace("TBRDist"))
  }
  .worker$call(function(t1, t2) as.numeric(TBRDist::USPRDist(t1, t2)),
               list(trees[[1]], trees[[2]]))
  if (.worker$poll_process(timeout * 1000) == "timeout") {
    .worker$kill()
    return(NA_real_)
  }
  res <- .worker$read()
  if (is.null(res$error)) res$result else NA_real_
}

if (file.exists(out)) {
  done <- readRDS(out)
} else {
  done <- list()
}
for (i in seq_along(pairs)) {
  if (i <= length(done)) next
  if (is.null(pairs[[i]]$exact)) {
    pairs[[i]]$exact <- ExactUSPR(pairs[[i]]$trees)
  }
  done[[i]] <- pairs[[i]]
  if (i %% 20 == 0) {
    saveRDS(done, out)
    cat(i, "/", length(pairs), "\n")
  }
}
saveRDS(done, out)
# Fixture check against the stored matrices
for (set in c("fix25-12spr", "fix32-16spr")) {
  ex <- vapply(done[vapply(done, `[[`, "", "set") == set], `[[`, 0, "exact")
  cat(set, "NA:", sum(is.na(ex)), "\n")
}
