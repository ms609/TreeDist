#' Generalized Robinson&ndash;Foulds distance
#' 
#' An internal function to calculate Generalized Robinson&ndash;Foulds
#' distances from splits.
#'
#' Note that no checks will be made to confirm that `splits1` and `splits2`
#' contain the same leaves in the same order.
#' This is the responsibility of the calling function.
#' 
#' @inheritParams SharedPhylogeneticInfoSplits
#' @param nTip Integer specifying the number of leaves in each split.
#' @param PairScorer function taking four arguments, `splits1`, `splits2`,
#' `nSplits1`, `nSplits2`, which should return the score of each pair of splits
#' in a two-dimensional matrix.  Additional parameters may be specified via 
#' \code{...}.
#' @param maximize Logical specifying whether the optimal matching maximizes
#' or minimizes the scores obtained by `PairScorer()`.
#' @param \dots Additional parameters to `PairScorer()`.
#' 
#' @return A numeric value specifying the score of the tree pairs under the 
#' specified pair scorer. If `reportMatching = TRUE`, attribute also list:
#' 
#' - `matching`: which split in `splits2` is optimally matched to each split in 
#'   `split1` (`NA` if not matched);
#' 
#' - `matchedSplits`: Textual representation of each match
#'
#' - `matchedScores`: Scores for matched split.
#' 
#' - `pairScores`: Calculated scores for each possible matching of each split.
#'
#' @keywords internal
#' @template MRS
#' @encoding UTF-8
#' @export
GeneralizedRF <- function(splits1, splits2, nTip, PairScorer, 
                           maximize, reportMatching, ...) {
  nSplits1 <- dim(splits1)[[1]]
  nSplits2 <- dim(splits2)[[1]]
  
  solution <- PairScorer(splits1, splits2, nTip, ...)
  ret <- solution[["score"]]
  
  if (reportMatching) {
    matching <- solution[["matching"]]
    matching[matching > nSplits2 | matching == 0L] <- NA
    if (nSplits1 < nSplits2) {
      matching <- matching[seq_len(nSplits1)]
    }
    attr(ret, "matching") <- matching
    
    # If reporting matching, we're not worried about performance
    pairScores <- matrix(0, nSplits1, nSplits2)
    for (i in seq_len(nSplits1)) {
      for (j in seq_len(nSplits2)) {
        pairScores[i, j] <- PairScorer(splits1[i, , drop = FALSE], 
                                       splits2[j, , drop = FALSE],
                                       nTip = nTip, ...)[["score"]]
      }
    }
    
    if (!is.null(attr(splits1, "tip.label"))) {
      matched1 <- !is.na(matching)
      matched2 <- matching[matched1]
      matched1 <- which(matched1)
      
      # Textual representation of matchings
      attr(ret, "matchedSplits") <- 
        ReportMatching(splits1[[matched1]],
                       splits2[[matched2]],
                       realMatch = if (maximize) {
                         pairScores[matrix(c(matched1, matched2), ncol = 2L)] > 0
                       } else rep(TRUE, length(matched1)))
    }
    
    attr(ret, "matchedScores") <- pairScores[
      matrix(c(seq_along(matching), matching), ncol = 2L)]
    
    attr(ret, "pairScores") <- pairScores
  }
  
  # Return:
  ret
}


#' @importFrom cli cli_progress_along
.MaxValue <- function(tree1, tree2, Value) {
  lab1 <- TipLabels(tree1)
  sameTips <- .AllTipsSame(lab1, TipLabels(tree2))
  if (!is.null(tree2)) {
    if (!sameTips) {
      trees <- .SharedOnly(tree1, tree2)
      tree1 <- trees[[1]]
      tree2 <- trees[[2]]
    }
  }
  
  if (is.null(tree2)) {
    if (sameTips) {
      value1 <- Value(tree1)
      # Much faster to discard unneeded than to only calculate required
      maxValue <- outer(value1, value1, "+")[, , drop = TRUE]
      maxValue[lower.tri(maxValue)]
    } else {
      # !sameTips implies that tree1 contains multiple trees
      pairs <- combn(seq_along(tree1), 2)
      nPairs <- dim(pairs)[[2]]
      
      vapply(cli_progress_along(seq_len(nPairs), "Calc max value"),
             function(pair) {
               i <- pairs[1, pair]
               j <- pairs[2, pair]
               common <- intersect(lab1[[i]], lab1[[j]])
               Value(KeepTip(tree1[[i]], common)) +
                 Value(KeepTip(tree1[[j]], common))
             }, double(1))
      
    }
  } else {
    value1 <- Value(tree1)
    if (sameTips) {
      outer(value1, Value(tree2), "+")[, , drop = TRUE]
    } else {
      value1 + Value(tree2)
    }
  }
}

.AllTipsSame <- function(x, y) {
  if (is.list(x)) {
    xPrime <- x[[1]]
    if (length(x) > 1 && !all(vapply(x[-1], setequal, logical(1), xPrime))) {
      return(FALSE)
    }
  } else {
    xPrime <- x
  }
  if (is.null(y)) {
    TRUE
  } else {
    if (is.list(y)) {
      all(vapply(y, setequal, logical(1), xPrime))
    } else {
      setequal(xPrime, y)
    }
  }
}

#' Are splits compatible?
#' 
#' Determine whether splits are compatible (concave); i.e. they can both occur
#' on a single tree.
#' 
#' @template split12Params
#' @return `SplitsCompatible()` returns a logical specifying whether the splits
#' provided are compatible with one another.
#' 
#' @examples 
#' A <- TRUE
#' B <- FALSE
#' SplitsCompatible(c(A, A, A, B, B, B),
#'                  c(A, A, B, B, B, B))
#' SplitsCompatible(c(A, A, A, B, B, B),
#'                  c(A, A, B, B, B, A))
#' @template MRS
#' @export
SplitsCompatible <- function(split1, split2) {
  # Return:
  (
    all (split1[split2]) ||
    all(split1[!split2]) ||
    all(!split1[split2]) ||
    all(!split1[!split2])
  )
}
