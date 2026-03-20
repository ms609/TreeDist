#' Consensus tree minimizing transfer distance
#'
#' Construct a consensus tree that minimizes the sum of transfer distances
#' to a set of input trees, using a greedy add-and-prune heuristic.
#'
#' Unlike the majority-rule consensus, which minimizes Robinson-Foulds
#' distance and can be highly unresolved when phylogenetic signal is low,
#' `TransferConsensus()` uses the finer-grained transfer distance
#' \insertCite{Lemoine2018}{TreeDist} to construct a more resolved consensus
#' tree.
#'
#' The algorithm pools all splits observed across input trees, computes
#' pairwise transfer distances between them, and greedily adds or removes
#' splits to minimize total transfer dissimilarity cost.  The approach
#' follows \insertCite{Takazawa2026;textual}{TreeDist}, reimplemented for
#' TreeDist's infrastructure.
#'
#' @param trees An object of class `multiPhylo`: the input trees.
#'   All trees must share the same tip labels.
#' @param scale Logical; if `TRUE` (default), use the scaled transfer
#'   distance (normalized by light-side size minus one).  If `FALSE`, use
#'   the unscaled (raw Hamming) transfer distance.
#' @param greedy Character string specifying the greedy strategy:
#'   `"best"` (default) picks the single highest-benefit action at each step;
#'   `"first"` picks the first improving action encountered (faster but
#'   potentially lower quality).
#' @param init Character string specifying the initial consensus:
#'   `"empty"` (default) starts with no splits (purely additive);
#'   `"majority"` starts with the majority-rule consensus and refines.
#'
#' @return A tree of class `phylo`.
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' library(TreeTools)
#' # Generate bootstrap-like trees
#' trees <- as.phylo(1:20, nTip = 12)
#'
#' # Transfer consensus (more resolved than majority rule)
#' tc <- TransferConsensus(trees)
#' plot(tc)
#'
#' # Compare resolution
#' mr <- Consensus(trees, p = 0.5)
#' cat("Majority-rule splits:", NSplits(mr), "\n")
#' cat("Transfer consensus splits:", NSplits(tc), "\n")
#'
#' @importFrom TreeTools as.Splits TipLabels NSplits Consensus StarTree
#' @export
TransferConsensus <- function(trees,
                              scale = TRUE,
                              greedy = c("best", "first"),
                              init = c("empty", "majority")) {
  greedy <- match.arg(greedy)
  init <- match.arg(init)

  if (!inherits(trees, "multiPhylo")) {
    stop("`trees` must be an object of class 'multiPhylo'.")
  }
  nTree <- length(trees)
  if (nTree < 2L) stop("Need at least 2 trees.")
  tipLabels <- TipLabels(trees[[1]])
  nTip <- length(tipLabels)

  # Convert each tree to a raw split matrix (TreeTools C++ internally)
  splitsList <- lapply(trees, function(tr) {
    sp <- as.Splits(tr, tipLabels)
    unclass(sp)
  })

  # Delegate all work to C++
  nThreads <- max(1L, getOption("TreeDist.threads",
                                 getOption("mc.cores", 1L)))
  res <- cpp_transfer_consensus(
    splitsList, nTip, scale,
    greedy_best_flag = (greedy == "best"),
    init_majority = (init == "majority"),
    n_threads = nThreads
  )

  included <- res$included
  if (!any(included)) {
    return(StarTree(tipLabels))
  }

  rawSplits <- res$raw_splits[included, , drop = FALSE]
  sp <- structure(rawSplits, nTip = nTip, tip.label = tipLabels,
                  class = "Splits")
  as.phylo(sp)
}


# ===========================================================================
# Internal helpers
# ===========================================================================

# popcount lookup for 0:255
.POPCOUNT <- as.integer(sapply(0:255, function(x) sum(as.integer(intToBits(x)))))

#' Pool unique splits, returning an integer (logical) matrix
#' @noRd
.PoolSplits <- function(trees, tipLabels) {
  nTip <- length(tipLabels)
  nTree <- length(trees)

  # Collect all splits as logical matrices and hash to unique set
  env <- new.env(hash = TRUE, parent = emptyenv())
  splitsList <- list()
  rawSplitsList <- list()
  countVec <- integer(0)
  treeMembers <- vector("list", nTree)
  nextIdx <- 1L

  for (i in seq_len(nTree)) {
    sp <- as.Splits(trees[[i]], tipLabels)
    logMat <- as.logical(sp)  # matrix: nSplits x nTip
    # The Splits object is a raw matrix internally; preserve structure
    rawMat <- unclass(sp)
    if (is.null(dim(logMat))) {
      logMat <- matrix(logMat, nrow = 1)
      rawMat <- matrix(rawMat, nrow = 1)
    }
    nSp <- nrow(logMat)
    # Canonicalise: ensure tip 1 is FALSE
    nSp <- nrow(logMat)
    members <- integer(nSp)
    for (j in seq_len(nSp)) {
      row <- logMat[j, ]
      rawRow <- rawMat[j, ]
      if (row[1]) {
        row <- !row
        rawRow <- .FlipRaw(rawRow, nTip)
      }
      key <- paste0(which(row), collapse = ",")
      idx <- env[[key]]
      if (is.null(idx)) {
        env[[key]] <- nextIdx
        splitsList[[nextIdx]] <- row
        rawSplitsList[[nextIdx]] <- rawRow
        countVec[nextIdx] <- 1L
        members[j] <- nextIdx
        nextIdx <- nextIdx + 1L
      } else {
        countVec[idx] <- countVec[idx] + 1L
        members[j] <- idx
      }
    }
    treeMembers[[i]] <- unique(members)
  }

  nSplits <- length(splitsList)
  splitMat <- matrix(FALSE, nrow = nSplits, ncol = nTip)
  nBytes <- length(rawSplitsList[[1]])
  rawMat <- matrix(as.raw(0), nrow = nSplits, ncol = nBytes)
  for (k in seq_len(nSplits)) {
    splitMat[k, ] <- splitsList[[k]]
    rawMat[k, ] <- rawSplitsList[[k]]
  }

  lightSide <- pmin(rowSums(splitMat), nTip - rowSums(splitMat))

  list(
    splits = splitMat,         # logical matrix nSplits x nTip
    rawSplits = rawMat,        # raw matrix nSplits x nBytes
    counts = countVec,
    lightSide = as.integer(lightSide),
    treeMembers = treeMembers
  )
}


.FlipRaw <- function(rawVec, nTip) {
  nBytes <- length(rawVec)
  usedBits <- ((nTip - 1L) %% 8L) + 1L
  lastMask <- as.raw(sum(2^(0:(usedBits - 1L))))
  out <- xor(rawVec, as.raw(0xff))
  out[nBytes] <- out[nBytes] & lastMask
  out
}


#' Pairwise transfer distance matrix using logical split matrix
#' Transfer distance = min(hamming(a, b), n - hamming(a, b))
#' hamming(a,b) when both are logical = sum(xor(a,b))
#' @noRd
.TransferDistMat <- function(splitMat, nTip) {
  # splitMat is logical: nSplits x nTip
  # hamming = number of differing positions
  # Use tcrossprod trick: hamming(a,b) = sum(a) + sum(b) - 2*sum(a&b)
  # = nAgreeing... actually let's just compute XOR directly.
  # Faster: hamming = nTip - 2 * (a %*% t(b) + (1-a) %*% t(1-b))
  #       = nTip - 2 * (a %*% t(b) + (nTip - rowSums(a) - colSums(b) + a %*% t(b)))
  # Simpler: agreement = a %*% t(b) + (1-a) %*% t(1-b)
  #                    = 2 * (a %*% t(b)) - rowSums(a) - rep(rowSums(b)) + nTip
  # hamming = nTip - agreement
  sm <- splitMat + 0L  # convert to integer
  ab <- tcrossprod(sm)  # sm %*% t(sm)
  rs <- rowSums(sm)
  hamming <- nTip - 2L * ab + outer(rs, rs, "+") - nTip
  # Simplifies to: hamming = outer(rs, rs, "+") - 2 * ab
  # Wait, let me re-derive:
  # agreement_ij = sum(a_i == b_j) = sum(a&b) + sum(!a & !b)
  #              = ab[i,j] + (nTip - rs[i] - rs[j] + ab[i,j])
  #              = 2*ab[i,j] + nTip - rs[i] - rs[j]
  # hamming_ij = nTip - agreement_ij = rs[i] + rs[j] - 2*ab[i,j]
  hamming <- outer(rs, rs, "+") - 2L * ab

  # Transfer distance = min(hamming, nTip - hamming)
  pmin(hamming, nTip - hamming)
}


#' Compute TD (transfer dissimilarity cost) for each split
#' @noRd
.ComputeTD <- function(DIST, sentDist, treeMembers, lightSide, nTree, scale) {
  nSplits <- nrow(DIST)
  TD <- numeric(nSplits)
  pMinus1 <- lightSide - 1L

  for (i in seq_len(nTree)) {
    idx <- treeMembers[[i]]
    # For each split b, min distance to any split in tree i
    if (length(idx) == 1L) {
      minD <- DIST[, idx]
    } else {
      minD <- apply(DIST[, idx, drop = FALSE], 1, min)
    }
    # Also consider sentinel distance
    minD <- pmin(minD, sentDist)
    if (scale) {
      TD <- TD + pmin(minD / pMinus1, 1)
    } else {
      TD <- TD + pmin(minD, pMinus1)
    }
  }
  TD
}


#' Compatibility matrix via vectorized bipartition check
#' @noRd
.CompatMat <- function(splitMat, nTip) {
  # Two splits compatible iff one of the 4 intersections (A&B, A&~B, ~A&B, ~A&~B) is empty
  sm <- splitMat + 0L  # integer 0/1 matrix
  notSm <- 1L - sm

  # A&B: a[i,] & b[j,] -> any nonzero -> tcrossprod(sm, sm) > 0
  ab   <- tcrossprod(sm, sm) > 0L
  anb  <- tcrossprod(sm, notSm) > 0L
  nab  <- tcrossprod(notSm, sm) > 0L
  nanb <- tcrossprod(notSm, notSm) > 0L

  # Compatible if at least one intersection is empty
  !ab | !anb | !nab | !nanb
}


#' Initialize MATCH/MATCH2 from currently included splits
#' @noRd
.InitMatches <- function(st, DIST, sentDist, lightSide, scale) {
  nSplits <- length(st$MATCH)
  incIdx <- which(st$incl)
  pMinus1 <- lightSide - 1L

  if (length(incIdx) == 0L) return(invisible())

  # For each split b, find 1st and 2nd closest among included
  subDIST <- DIST[, incIdx, drop = FALSE]
  for (b in seq_len(nSplits)) {
    dists <- subDIST[b, ]
    ord <- order(dists)
    bestDist <- dists[ord[1]]
    thresh <- pMinus1[b]
    if (scale && bestDist / thresh >= 1) next
    if (!scale && bestDist >= thresh) next
    st$MATCH[b] <- incIdx[ord[1]]
    if (length(ord) > 1L) {
      secDist <- dists[ord[2]]
      if (scale && secDist / thresh < 1) {
        st$MATCH2[b] <- incIdx[ord[2]]
      } else if (!scale && secDist < thresh) {
        st$MATCH2[b] <- incIdx[ord[2]]
      }
    }
  }
}


#' Get distance from split b to its match (NA = sentinel)
#' @noRd
.Dist <- function(b, idx, DIST, sentDist) {
  if (is.na(idx)) sentDist[b] else DIST[b, idx]
}


#' Vectorized add-benefit: returns benefit for each candidate
#' @noRd
.AddBenefitVec <- function(candidates, st, DIST, sentDist, TD, counts,
                           lightSide, scale) {
  nSplits <- length(st$MATCH)
  pMinus1 <- lightSide - 1L

  # Current match distances for all splits
  matchDist <- ifelse(is.na(st$MATCH), sentDist, .DiagDist(st$MATCH, DIST, sentDist))

  benefits <- numeric(length(candidates))
  for (ci in seq_along(candidates)) {
    c <- candidates[ci]
    newDist <- DIST[, c]
    if (scale) {
      diff <- (matchDist - newDist) / pMinus1
    } else {
      diff <- matchDist - newDist
    }
    diff[diff < 0] <- 0
    benefits[ci] <- sum(diff * counts) - TD[c]
  }
  benefits
}

#' Helper: get DIST\[b, MATCH\[b\]\] for all b, vectorized
#' @noRd
.DiagDist <- function(matchVec, DIST, sentDist) {
  n <- length(matchVec)
  out <- numeric(n)
  notNA <- !is.na(matchVec)
  if (any(notNA)) {
    out[notNA] <- DIST[cbind(which(notNA), matchVec[notNA])]
  }
  out[!notNA] <- sentDist[!notNA]
  out
}


#' Vectorized remove-benefit
#' @noRd
.RemoveBenefitVec <- function(candidates, st, DIST, sentDist, TD, counts,
                              lightSide, scale) {
  nSplits <- length(st$MATCH)
  pMinus1 <- lightSide - 1L

  # For remove, only splits whose MATCH == candidate are affected
  benefits <- numeric(length(candidates))
  matchDist <- .DiagDist(st$MATCH, DIST, sentDist)
  match2Dist <- .DiagDist(st$MATCH2, DIST, sentDist)

  for (ci in seq_along(candidates)) {
    c <- candidates[ci]
    affected <- st$MATCH == c & !is.na(st$MATCH)
    if (any(affected)) {
      if (scale) {
        fn_cost <- sum((DIST[affected, c] - match2Dist[affected]) /
                         pMinus1[affected] * counts[affected])
      } else {
        fn_cost <- sum((DIST[affected, c] - match2Dist[affected]) *
                         counts[affected])
      }
    } else {
      fn_cost <- 0
    }
    benefits[ci] <- TD[c] + fn_cost
  }
  benefits
}


.DoAdd <- function(branchIdx, st, DIST, sentDist) {
  st$incl[branchIdx] <- TRUE
  nSplits <- length(st$MATCH)

  curMatchDist <- .DiagDist(st$MATCH, DIST, sentDist)
  newDist <- DIST[, branchIdx]

  better <- newDist < curMatchDist
  if (any(better)) {
    st$MATCH2[better] <- st$MATCH[better]
    st$MATCH[better] <- branchIdx
  }

  # Check if it becomes second match for others
  notBetter <- !better
  secDist <- .DiagDist(st$MATCH2, DIST, sentDist)
  betterSec <- notBetter & (newDist < secDist)
  if (any(betterSec)) {
    st$MATCH2[betterSec] <- branchIdx
  }
}


.DoRemove <- function(branchIdx, st, DIST, sentDist, lightSide, scale) {
  st$incl[branchIdx] <- FALSE
  nSplits <- length(st$MATCH)
  curInc <- which(st$incl)
  pMinus1 <- lightSide - 1L

  # Splits whose first match was branchIdx
  affected1 <- which(st$MATCH == branchIdx & !is.na(st$MATCH))
  if (length(affected1)) {
    st$MATCH[affected1] <- st$MATCH2[affected1]
    for (b in affected1) {
      if (is.na(st$MATCH[b])) {
        # Promoted value was sentinel — rescan for actual closest
        st$MATCH[b] <- .FindSecond(b, NA_integer_, curInc, DIST,
                                    pMinus1, scale)
      }
      # Find new second match
      st$MATCH2[b] <- .FindSecond(b, st$MATCH[b], curInc, DIST,
                                  pMinus1, scale)
    }
  }

  # Splits whose second match was branchIdx
  affected2 <- which(st$MATCH2 == branchIdx & !is.na(st$MATCH2))
  if (length(affected2)) {
    for (b in affected2) {
      st$MATCH2[b] <- .FindSecond(b, st$MATCH[b], curInc, DIST,
                                  pMinus1, scale)
    }
  }
}


.FindSecond <- function(b, matchIdx, curInc, DIST, pMinus1, scale) {
  cands <- if (is.na(matchIdx)) curInc else setdiff(curInc, matchIdx)
  if (length(cands) == 0L) return(NA_integer_)
  dists <- DIST[b, cands]
  best <- cands[which.min(dists)]
  bestDist <- DIST[b, best]
  if (scale && bestDist / pMinus1[b] >= 1) return(NA_integer_)
  if (!scale && bestDist >= pMinus1[b]) return(NA_integer_)
  best
}


.IsCompat <- function(idx, incl, compat, nTip) {
  incIdx <- which(incl)
  if (length(incIdx) == 0L) return(TRUE)
  if (length(incIdx) >= nTip - 3L) return(FALSE)
  all(compat[idx, incIdx])
}


.GreedyBest <- function(st, DIST, sentDist, TD, counts, lightSide,
                        compat, sortOrd, scale, nSplits, nTip) {
  repeat {
    # Evaluate all candidates
    addCands <- integer(0)
    remCands <- integer(0)
    for (idx in sortOrd) {
      if (st$incl[idx]) {
        remCands <- c(remCands, idx)
      } else if (.IsCompat(idx, st$incl, compat, nTip)) {
        addCands <- c(addCands, idx)
      }
    }

    bestBen <- 0
    bestIdx <- 0L
    bestAction <- ""

    if (length(addCands)) {
      addBen <- .AddBenefitVec(addCands, st, DIST, sentDist, TD, counts,
                               lightSide, scale)
      mx <- max(addBen)
      if (mx > bestBen) {
        bestBen <- mx
        bestIdx <- addCands[which.max(addBen)]
        bestAction <- "add"
      }
    }
    if (length(remCands)) {
      remBen <- .RemoveBenefitVec(remCands, st, DIST, sentDist, TD, counts,
                                  lightSide, scale)
      mx <- max(remBen)
      if (mx > bestBen) {
        bestBen <- mx
        bestIdx <- remCands[which.max(remBen)]
        bestAction <- "remove"
      }
    }

    if (bestBen <= 0) break
    if (bestAction == "add") {
      .DoAdd(bestIdx, st, DIST, sentDist)
    } else {
      .DoRemove(bestIdx, st, DIST, sentDist, lightSide, scale)
    }
  }
}


.GreedyFirst <- function(st, DIST, sentDist, TD, counts, lightSide,
                         compat, sortOrd, scale, nSplits, nTip) {
  improving <- TRUE
  while (improving) {
    improving <- FALSE
    matchDist <- .DiagDist(st$MATCH, DIST, sentDist)
    match2Dist <- .DiagDist(st$MATCH2, DIST, sentDist)
    pMinus1 <- lightSide - 1L

    for (idx in sortOrd) {
      if (st$incl[idx]) {
        # Quick remove benefit
        affected <- st$MATCH == idx & !is.na(st$MATCH)
        if (any(affected)) {
          if (scale) {
            fn <- sum((DIST[affected, idx] - match2Dist[affected]) /
                        pMinus1[affected] * counts[affected])
          } else {
            fn <- sum((DIST[affected, idx] - match2Dist[affected]) *
                        counts[affected])
          }
        } else {
          fn <- 0
        }
        if (TD[idx] + fn > 0) {
          .DoRemove(idx, st, DIST, sentDist, lightSide, scale)
          improving <- TRUE
          break
        }
      } else if (.IsCompat(idx, st$incl, compat, nTip)) {
        newDist <- DIST[, idx]
        if (scale) {
          diff <- pmax((matchDist - newDist) / pMinus1, 0)
        } else {
          diff <- pmax(matchDist - newDist, 0)
        }
        if (sum(diff * counts) - TD[idx] > 0) {
          .DoAdd(idx, st, DIST, sentDist)
          improving <- TRUE
          break
        }
      }
    }
  }
}


.SplitsToPhylo <- function(rawSplits, included, tipLabels, nTip) {
  selectedIdx <- which(included)
  if (length(selectedIdx) == 0L) {
    return(TreeTools::StarTree(tipLabels))
  }
  selectedRaw <- rawSplits[selectedIdx, , drop = FALSE]
  sp <- structure(selectedRaw, nTip = nTip, tip.label = tipLabels,
                  class = "Splits")
  as.phylo(sp)
}
