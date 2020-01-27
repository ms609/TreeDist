#' Generalized Robinson-Foulds distance
#' 
#' An internal function to calculate Generalized Robinson-Foulds distance from
#' splits.
#'
#' Note that no checks will be made to confirm that splits1 and splits2 contain
#' the same leaves in the same order.  This is the responsibility of the calling
#' function.
#' 
#' @inheritParams SharedPhylogeneticInfoSplits
#' @template nTipParam
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
#' @export
#' @references \insertRef{Jonker1987}{TreeDist}
GeneralizedRF <- function (splits1, splits2, nTip, PairScorer, 
                           maximize, reportMatching, ...) {
  nSplits1 <- dim(splits1)[1]
  nSplits2 <- dim(splits2)[1]
  if (nSplits1 == 0 || nSplits2 == 0) return (0L)
  
  solution <- PairScorer(splits1, splits2, nTip, ...)
  ret <- solution$score
  
  if (reportMatching) {
    matching <- solution$matching
    matching[matching > nSplits2] <- NA
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
                         pairScores[matrix(c(matched1, matched2), ncol=2L)] > 0
                       } else TRUE)
    }
  }
  # Return:
  ret
}

#' Are splits compatible?
#' 
#' Splits are compatible if they are concave; i.e. they can both be true
#' simultaneously.
#' 
#' @template split12Params
#' @return `SplitsCompatible` returns a logical specifying whether the splits
#' provided are compatible with one another.
#' 
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
