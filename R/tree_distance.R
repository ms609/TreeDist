#' Generalized Robinson-Foulds distance
#' 
#' An internal function to calculate Generalized Robinson-Foulds distances from
#' splits.
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
#' \dots.
#' @param \dots Additional parameters to `PairScorer`
#' 
#' @return A numeric value specifying the score of the tree pairs under the 
#' specified pair scorer. If `reportMatching = TRUE`, attribute also list:
#' 
#' - `matching`: which split in `splits2` is optimally matched to each split in 
#' `split1` (`NA` if not matched);
#'
#' - `pairScores`: Calculated scores for each possible matching of each split.
#' 
#' - `matchedSplits`: Textual representation of each match
#' 
#' @keywords internal
#' @template MRS
#' @encoding UTF-8
#' @export
#' @references \insertRef{Jonker1987}{TreeDist}
GeneralizedRF <- function (splits1, splits2, nTip, PairScorer, 
                           maximize, reportMatching, ...) {
  nSplits1 <- dim(splits1)[1]
  nSplits2 <- dim(splits2)[1]
  
  solution <- PairScorer(splits1, splits2, nTip, ...)
  ret <- solution$score
  
  if (reportMatching) {
    matching <- solution$matching
    matching[matching > nSplits2 | matching == 0L] <- NA
    if (nSplits1 < nSplits2) {
      matching <- matching[seq_len(nSplits1)]
    }
    attr(ret, 'matching') <- matching
    
    # If reporting matching, we're not worried about performance
    pairScores <- matrix(0, nSplits1, nSplits2)
    for (i in seq_len(nSplits1)) {
      for (j in seq_len(nSplits2)) {
        pairScores[i, j] <- PairScorer(splits1[i, , drop = FALSE], 
                                       splits2[j, , drop = FALSE],
                                       nTip = nTip, ...)$score
      }
    }
    attr(ret, 'pairScores') <- pairScores
    
    if (!is.null(attr(splits1, 'tip.label'))) {
      matched1 <- !is.na(matching)
      matched2 <- matching[matched1]
      matched1 <- which(matched1)
      
      attr(ret, 'matchedSplits') <- 
        ReportMatching(splits1[[matched1]],
                       splits2[[matched2]],
                       realMatch = if (maximize) {
                         pairScores[matrix(c(matched1, matched2), ncol = 2L)] > 0
                       } else rep(TRUE, length(matched1)))
    }
  }
  # Return:
  ret
}

.MaxValue <- function (tree1, tree2, Value) {
  value1 <- Value(tree1)
  if (is.null(tree2)) {
    maxValue <- outer(value1, value1, '+')[, , drop = TRUE]
    maxValue[lower.tri(maxValue)]
  } else {
    outer(value1, Value(tree2), '+')[, , drop = TRUE]
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
SplitsCompatible <- function (split1, split2) {
  # Return:
  (
    all (split1[split2]) ||
    all(split1[!split2]) ||
    all(!split1[split2]) ||
    all(!split1[!split2])
  )
}
