# Leaf selectors for SPRDist(method = "Rogue"), with a harness that copies the
# depth-1 loop of .SPRRogue() but takes the choice of leaf from `scorer`.
#
# Requires devtools::load_all() of the TreeDist worktree (internal functions).

suppressMessages({
  library("TreeTools")
})

# Log of each scoring call: number of leaves, size of the Fitch tie set, etc.
.log <- new.env()
ResetLog <- function() .log$rows <- list()
LogRow <- function(...) .log$rows[[length(.log$rows) + 1L]] <- list(...)
ResetLog()

# Copy of .SPRRogue() (depth 1, C++ path), with the leaf choice delegated.
# scorer(tr1, tr2) returns the index of the leaf to drop.
RogueWith <- function(tree1, tree2, scorer, check = TRUE) {
  moves <- 0
  reduced <- ReduceTrees(tree1, tree2, check = check)
  while (!is.null(reduced)) {
    tr1 <- reduced[[1]]
    tr2 <- reduced[[2]]
    nTip <- NTip(tr1)
    if (nTip == 4) {
      return(moves + 1)
    }
    if (nTip < 8) {
      sp1 <- as.Splits(tr1)
      sp2 <- as.Splits(tr2, tr1)
      return(moves + switch(nTip - 4L,
                            .SPRExact5(sp1, sp2),
                            .SPRExact6(sp1, sp2),
                            .SPRExact7(sp1, sp2)))
    }
    subtips1 <- spr_first_shared_split(tr1, tr2)
    if (length(subtips1)) {
      subtips2 <- !subtips1
      subtips1[!subtips1][[1]] <- TRUE
      subtips2[!subtips2][[1]] <- TRUE
      return(moves +
               RogueWith(KeepTipPostorder(tr1, subtips1),
                         KeepTipPostorder(tr2, subtips1), scorer) +
               RogueWith(KeepTipPostorder(tr1, subtips2),
                         KeepTipPostorder(tr2, subtips2), scorer))
    }
    drop <- logical(nTip)
    drop[[scorer(tr1, tr2)]] <- TRUE
    reduced <- spr_keep_and_reduce(tr1, tr2, !drop)
    if (is.null(reduced[[1]])) {
      reduced <- NULL
    }
    moves <- moves + 1
  }
  moves
}

# Current selector: least bounded-MCI distance after dropping each leaf.
MCIScores <- function(tr1, tr2) spr_rogue_scores(tr1, tr2, 1L)[["scores"]]
ScorerMCI <- function(tr1, tr2) {
  s <- MCIScores(tr1, tr2)
  LogRow(n = NTip(tr1), sel = "mci")
  which.min(s)
}

# Incremental Fitch for binary characters.
#
# For each leaf x, returns the summed length of `tree` over the characters in
# `chars` (logical, nChar x nTip) when x is coded `?` in every character.
# Coding a leaf `?` is equivalent to deleting it, so length(T, x = ?) is the
# Fitch cost of T minus x's subtree, which one uppass gives for every x:
#   D[v], dc[v]: downpass state set and cost of the subtree at v;
#   U[v], uc[v]: set and cost of T minus subtree(v), seen from v's attachment.
#   Root children a, b: U[a] = D[b], uc[a] = dc[b].
#   v with parent p, sibling s: U[v] = F(D[s], U[p]),
#     uc[v] = dc[s] + uc[p] + [D[s] & U[p] empty].
# State sets are bitmasks (1 = {0}, 2 = {1}), vectorised over characters.
# Cost: O(nNode) vector ops of length nChar, i.e. O(n * nChar).
FitchAmbiguous <- function(tree, chars) {
  tree <- Preorder(RootTree(tree, 1))
  edge <- tree[["edge"]]
  nTip <- NTip(tree)
  nNode <- nTip + tree[["Nnode"]]
  root <- nTip + 1L
  parent <- integer(nNode)
  parent[edge[, 2]] <- edge[, 1]
  kids <- split(edge[, 2], edge[, 1])
  nChar <- dim(chars)[[1]]
  D <- matrix(0L, nChar, nNode)
  dc <- matrix(0L, nChar, nNode)
  D[, seq_len(nTip)] <- 1L + chars
  # Postorder: reverse of preorder edge list
  for (v in rev(unique(edge[, 1]))) {
    k <- kids[[as.character(v)]]
    stopifnot(length(k) == 2L)
    a <- k[[1]]; b <- k[[2]]
    I <- bitwAnd(D[, a], D[, b])
    empty <- I == 0L
    D[, v] <- ifelse(empty, bitwOr(D[, a], D[, b]), I)
    dc[, v] <- dc[, a] + dc[, b] + empty
  }
  U <- matrix(0L, nChar, nNode)
  uc <- matrix(0L, nChar, nNode)
  for (v in edge[, 2]) { # Preorder: parents before children
    p <- parent[[v]]
    k <- kids[[as.character(p)]]
    s <- k[[which(k != v)]]
    if (p == root) {
      U[, v] <- D[, s]
      uc[, v] <- dc[, s]
    } else {
      I <- bitwAnd(D[, s], U[, p])
      empty <- I == 0L
      U[, v] <- ifelse(empty, bitwOr(D[, s], U[, p]), I)
      uc[, v] <- dc[, s] + uc[, p] + empty
    }
  }
  list(length = sum(dc[, root]),
       ambiguous = colSums(uc[, seq_len(nTip), drop = FALSE]))
}

# Splits of `tr2` as binary characters ordered by tr1's tips (nChar x nTip).
SplitChars <- function(trA, trB) {
  as.logical(as.Splits(trB, tipLabels = TipLabels(trA)))
}

# Length added by each leaf: steps saved when it is coded `?` everywhere.
# direction: "both", "12" (tree2's splits on tree1) or "21".
FitchScores <- function(tr1, tr2, direction = "both") {
  score <- numeric(NTip(tr1))
  if (direction %in% c("both", "12")) {
    f <- FitchAmbiguous(tr1, SplitChars(tr1, tr2))
    score <- score + f$length - f$ambiguous
  }
  if (direction %in% c("both", "21")) {
    # tr2's tip order matches tr1's after ReduceTrees(); index by label anyway.
    f <- FitchAmbiguous(tr2, SplitChars(tr2, tr1))
    add <- f$length - f$ambiguous
    score <- score + add[match(TipLabels(tr1), TipLabels(tr2))]
  }
  score
}

# Pure Fitch: drop the leaf adding most length; ties by leaf order.
# k = 0: no MCI; k > 0: MCI among the k highest-scoring leaves (ties by order);
# k = "tied": MCI among all leaves tied at the maximum.
# last = TRUE breaks Fitch ties by the last leaf instead (order sensitivity).
MakeFitchScorer <- function(direction = "both", k = 0, last = FALSE) {
  force(direction); force(k); force(last)
  function(tr1, tr2) {
    fs <- FitchScores(tr1, tr2, direction)
    nTied <- sum(fs == max(fs))
    ord <- if (last) rev(order(fs)) else order(-fs)
    cand <- if (identical(k, "tied")) {
      which(fs == max(fs))
    } else if (k > 0) {
      ord[seq_len(min(k, length(fs)))]
    } else {
      ord[[1]]
    }
    LogRow(n = NTip(tr1), sel = paste(direction, k), nTied = nTied,
           nCand = length(cand))
    if (length(cand) == 1L) {
      cand
    } else {
      mci <- MCIScores(tr1, tr2)
      cand[[which.min(mci[cand])]]
    }
  }
}

# Fitch score as a tie-breaker *after* MCI (for comparison only)
ScorerMCIThenFitch <- function(tr1, tr2) {
  mci <- MCIScores(tr1, tr2)
  best <- which(mci == min(mci))
  if (length(best) == 1L) return(best)
  fs <- FitchScores(tr1, tr2, "both")
  best[[which.max(fs[best])]]
}

# Full MCI scoring on reduced trees of <= nMax leaves; Fitch above.
MakeHybridScorer <- function(nMax = 20, k = 0) {
  fitch <- MakeFitchScorer("both", k)
  function(tr1, tr2) {
    if (NTip(tr1) <= nMax) ScorerMCI(tr1, tr2) else fitch(tr1, tr2)
  }
}
