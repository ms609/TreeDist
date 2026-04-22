#' Transfer dissimilarity between phylogenetic trees
#'
#' Compute the transfer dissimilarity between phylogenetic trees, as defined
#' by \insertCite{Takazawa2026;textual}{TreeDist}.
#' The transfer dissimilarity uses the transfer distance
#' \insertCite{Lemoine2018}{TreeDist} to compare bipartitions, providing a
#' finer-grained measure than the Robinson--Foulds distance.  Each branch in
#' each tree is scored by how many taxa must be moved to match its closest
#' counterpart in the other tree, and these scores are summed.
#'
#' The `scaled` variant divides each branch's contribution by its depth minus
#' one, giving equal weight to all branches regardless of their depth (analogous
#' to the Robinson--Foulds distance).  The `unscaled` variant uses raw transfer
#' distances, giving more weight to deep branches.  Neither variant satisfies
#' the triangle inequality for trees with six or more tips.
#'
#' @inheritParams TreeDistance
#' @param scale Logical; if `TRUE` (default), use the scaled transfer
#'   dissimilarity.  If `FALSE`, use the unscaled transfer dissimilarity.
#'
#' @return `TransferDist()` returns an object of class `dist` (if `tree2` is
#'   `NULL`), a numeric matrix (if both `tree1` and `tree2` are lists), or a
#'   numeric value (for a single pair).  If `reportMatching = TRUE`, the
#'   return value carries `matching` and `pairScores` attributes.
#'
#' @section Normalization:
#'
#' When `normalize = TRUE`, the scaled transfer dissimilarity is divided by
#' `2 * (n - 3)`, placing it in the range \[0, 1\].  The unscaled version is
#' divided by the maximum possible unscaled dissimilarity
#' (following \insertCite{Takazawa2026;textual}{TreeDist}).
#'
#' @examples
#' library(TreeTools)
#' TransferDist(BalancedTree(8), PectinateTree(8))
#' TransferDist(BalancedTree(8), PectinateTree(8), scale = FALSE)
#'
#' # All-pairs
#' TransferDist(as.phylo(0:5, 8))
#'
#' @references
#' \insertAllCited{}
#'
#' @family tree distances
#'
#' @importFrom TreeTools as.Splits TipLabels NSplits
#' @export
TransferDist <- function(tree1, tree2 = NULL, scale = TRUE,
                         normalize = FALSE, reportMatching = FALSE) {
  
  # --- Fast path: all-pairs (tree2 = NULL) ---
  if (is.null(tree2) && !reportMatching) {
    fast <- .TransferDistAllPairs(tree1, scale)
    if (!is.null(fast)) {
      if (!isFALSE(normalize)) {
        nTip <- length(TipLabels(tree1[[1]]))
        denom <- .TransferNormDenom(nTip, scale)
        fast <- fast / denom
      }
      return(fast)
    }
  }
  
  # --- Fast path: cross-pairs ---
  if (!is.null(tree2) && !reportMatching) {
    fast <- .TransferDistCrossPairs(tree1, tree2, scale)
    if (!is.null(fast)) {
      if (!isFALSE(normalize)) {
        nTip <- length(TipLabels(
          if (inherits(tree1, c("phylo", "Splits"))) tree1 else tree1[[1]]))
        denom <- .TransferNormDenom(nTip, scale)
        fast <- fast / denom
      }
      return(fast)
    }
  }
  
  # --- Generic fallback via CalculateTreeDistance ---
  # Capture `scale` in the closure for the Splits-level function
  Func <- function(splits1, splits2, nTip, reportMatching = FALSE) {
    TransferDistSplits(splits1, splits2, nTip = nTip,
                       scale = scale, reportMatching = reportMatching)
  }
  unnormalized <- CalculateTreeDistance(Func, tree1, tree2, reportMatching)
  
  if (!isFALSE(normalize)) {
    nTip <- length(TipLabels(
      if (inherits(tree1, c("phylo", "Splits"))) tree1 else tree1[[1]]))
    denom <- .TransferNormDenom(nTip, scale)
    unnormalized <- unnormalized / denom
  }
  
  unnormalized
}

#' @rdname TransferDist
#' @export
TransferDistance <- TransferDist

#' @rdname TransferDist
#' @inheritParams SharedPhylogeneticInfoSplits
#' @export
TransferDistSplits <- function(splits1, splits2,
                                nTip = attr(splits1, "nTip"),
                                scale = TRUE,
                                reportMatching = FALSE) {
  .ValidateSplitArgs(splits1, splits2, nTip)
  solution <- cpp_transfer_dist_scored(splits1, splits2,
                                        nTip = as.integer(nTip),
                                        scale = scale)
  ret <- solution[["score"]]
  
  if (reportMatching) {
    nSplits1 <- nrow(splits1)
    nSplits2 <- nrow(splits2)
    matching <- solution[["matching"]]
    matching[matching > nSplits2 | matching == 0L] <- NA
    if (nSplits1 < nSplits2) {
      matching <- matching[seq_len(nSplits1)]
    }
    attr(ret, "matching") <- matching
    
    # Compute full pairwise score matrix for reportMatching
    pairScores <- matrix(0, nSplits1, nSplits2)
    for (i in seq_len(nSplits1)) {
      for (j in seq_len(nSplits2)) {
        # Per-pair: the transfer distance contribution
        # This is the individual δ(b_i, b*_j) / (depth(b_i) - 1) for scaled
        # or min(δ(b_i, b*_j), depth(b_i) - 1) for unscaled
        pair_res <- cpp_transfer_dist_scored(
          splits1[i, , drop = FALSE],
          splits2[j, , drop = FALSE],
          nTip = as.integer(nTip),
          scale = scale
        )
        pairScores[i, j] <- pair_res[["score"]]
      }
    }
    
    if (!is.null(attr(splits1, "tip.label"))) {
      matched1 <- !is.na(matching)
      matched2 <- matching[matched1]
      matched1 <- which(matched1)
      attr(ret, "matchedSplits") <-
        ReportMatching(splits1[[matched1]], splits2[[matched2]],
                       realMatch = rep(TRUE, length(matched1)))
    }
    
    attr(ret, "matchedScores") <- pairScores[
      matrix(c(seq_along(matching), matching), ncol = 2L)]
    attr(ret, "pairScores") <- pairScores
  }
  
  ret
}


# ============================================================================
# Internal helpers
# ============================================================================

# All-pairs fast path
.TransferDistAllPairs <- function(tree1, scale) {
  if (inherits(tree1, c("phylo", "Splits"))) return(NULL)
  if (length(tree1) < 2L) return(NULL)
  
  tipLabels <- TipLabels(tree1[[1]])
  if (is.null(tipLabels)) return(NULL)
  nTip <- length(tipLabels)
  if (nTip < 4L) return(NULL)
  .CheckMaxTips(nTip)
  
  # Check all trees share same tip set
  allLabels <- TipLabels(tree1)
  if (is.list(allLabels)) {
    if (!all(vapply(allLabels[-1], setequal, logical(1), tipLabels))) {
      return(NULL)
    }
  }
  
  splitsList <- lapply(tree1, function(tr) {
    unclass(as.Splits(tr, tipLabels))
  })
  
  nThreads <- max(1L, getOption("TreeDist.threads",
                                 getOption("mc.cores", 1L)))
  
  dists <- cpp_transfer_dist_all_pairs(splitsList, nTip, scale, nThreads)
  
  N <- length(tree1)
  attr(dists, "Size") <- N
  attr(dists, "Diag") <- FALSE
  attr(dists, "Upper") <- FALSE
  nms <- names(tree1)
  if (!is.null(nms)) attr(dists, "Labels") <- nms
  class(dists) <- "dist"
  dists
}


# Cross-pairs fast path
.TransferDistCrossPairs <- function(tree1, tree2, scale) {
  single1 <- inherits(tree1, c("phylo", "Splits"))
  single2 <- inherits(tree2, c("phylo", "Splits"))
  if (single1 && single2) return(NULL) # use generic path for single pair
  
  trees1 <- if (single1) list(tree1) else tree1
  trees2 <- if (single2) list(tree2) else tree2
  
  tipLabels <- TipLabels(trees1[[1]])
  if (is.null(tipLabels)) return(NULL)
  nTip <- length(tipLabels)
  if (nTip < 4L) return(NULL)
  .CheckMaxTips(nTip)
  
  # Check all trees share same tip set
  allLabels1 <- TipLabels(trees1)
  allLabels2 <- TipLabels(trees2)
  if (is.list(allLabels1)) {
    if (!all(vapply(allLabels1, setequal, logical(1), tipLabels))) return(NULL)
  }
  if (is.list(allLabels2)) {
    if (!all(vapply(allLabels2, setequal, logical(1), tipLabels))) return(NULL)
  } else {
    if (!setequal(allLabels2, tipLabels)) return(NULL)
  }
  
  splits1 <- lapply(trees1, function(tr) unclass(as.Splits(tr, tipLabels)))
  splits2 <- lapply(trees2, function(tr) unclass(as.Splits(tr, tipLabels)))
  
  nThreads <- max(1L, getOption("TreeDist.threads",
                                 getOption("mc.cores", 1L)))
  
  mat <- cpp_transfer_dist_cross_pairs(splits1, splits2, nTip, scale, nThreads)
  
  rownames(mat) <- names(trees1)
  colnames(mat) <- names(trees2)
  
  # If one input was a single tree, simplify to vector
  if (single1) return(mat[1, ])
  if (single2) return(mat[, 1])
  mat
}


# Normalization denominator
.TransferNormDenom <- function(nTip, scale) {
  if (scale) {
    # Scaled: each tree contributes at most (n-3) branches × 1.0
    2.0 * (nTip - 3L)
  } else {
    # Unscaled: maximum dissimilarity between two caterpillar trees
    # Takazawa 2026: (n^2 - 2n + 4)/4 for even n, (n^2 - 2n + 5)/4 for odd n
    if (nTip %% 2L == 0L) {
      (nTip^2 - 2 * nTip + 4) / 4
    } else {
      (nTip^2 - 2 * nTip + 5) / 4
    }
  }
}
