#' Generalized Robinson-Foulds distance
#' 
#' Calculate Generalized Robinson-Foulds distance from splits.
#' 
#' @inheritParams MutualPhylogeneticInfoSplits
#' @param PairScorer function taking four arguments, `splits1`, `splits2`,
#' `nSplits1`, `nSplits2`, which should return the score of each pair of splits
#' in a two-dimensional matrix.  Additional parameters may be specified via 
#' \dots.
#' @param \dots Additional parameters to `PairScorer`
#' 
#' @return The results of `TreeDistanceReturn` under the parameters provided
#' 
#' @keywords internal
#' @author Martin R. Smith
#' @export
GeneralizedRF <- function (splits1, splits2, PairScorer, 
                           maximize, reportMatching, ...) {
  dimSplits1 <- dim(splits1)
  dimSplits2 <- dim(splits2)
  nSplits1 <- dimSplits1[2]
  nSplits2 <- dimSplits2[2]
  if (nSplits1 == 0 || nSplits2 == 0) return (0L)
  nTerminals <- dimSplits1[1]
  if (dimSplits2[1] != nTerminals) {
    stop("Split rows must bear identical labels")
  }
  
  taxonNames1 <- rownames(splits1)
  taxonNames2 <- rownames(splits2)
  
  if (!is.null(taxonNames2)) {
    splits2 <- unname(splits2[taxonNames1, , drop=FALSE])
    splits1 <- unname(splits1) # split2[split1] faster without names
  }
  
  pairScores <- PairScorer(splits1, splits2, nSplits1, nSplits2, ...)
  
  # Return:
  TreeDistanceReturn(pairScores, maximize, reportMatching, 
                     splits1, splits2, taxonNames1)
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
#' @author Martin R. Smith
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

#' Tree Distance Return
#' 
#' Generates the return value for Generalized Robinson-Foulds distances,
#' with (optionally) a report of the matching.
#' 
#' @param pairScores Two-dimensional array listing the score of each possible
#' pairing of `splits1` (rows) with `splits2` (columns).
#' @param reportMatching Logical specifying whether to report details
#' of an optimal matching.
#' @template splits12params
#' @param taxonNames Character vector listing the names corresponding to each
#' row of `splits1` and `splits2`.
#' Corresponding parameters from within each function
#' 
#' @keywords internal
#' @author Martin R. Smith
#' @importFrom clue solve_LSAP
#' @export
TreeDistanceReturn <- function (pairScores, maximize = FALSE,
                                reportMatching,  
                                splits1, splits2,
                                taxonNames = NULL) {
  dimScores <- dim(pairScores)
  
  if (dimScores[1] > dimScores[2]) {
    if (dimScores[2] == 1) {
      Optimal <- if (maximize) which.max else which.min
      solution <- Optimal(pairScores)
    } else {
      solution <- solve_LSAP(t(pairScores), maximize)
    }
    optimalMatching <- structure(match(seq_len(dimScores[1]), solution),
                                 class='solve_LSAP')
  } else {
    if (dimScores[1] == 1) {
      Optimal <- if (maximize) which.max else which.min
      optimalMatching <- structure(Optimal(pairScores), class='solve_LSAP')
    } else {
      optimalMatching <- solve_LSAP(pairScores, maximize)
    }
  }
  
  
  matched <- !is.na(optimalMatching)
  matched1 <- which(matched)
  matched2 <- optimalMatching[matched]
  
  ret <- sum(pairScores[matrix(c(matched1, matched2), ncol=2L)])
  
  if (reportMatching) {
    attributes(ret) <- list(
      pairScores = pairScores,
      matching = optimalMatching
    )
    
    if (!is.null(taxonNames)) {
      attr(ret, 'matchedSplits') <- 
        ReportMatching(splits1[, matched1, drop = FALSE], 
                       splits2[, matched2, drop = FALSE],
                       taxonNames)
    }
  }
  # Return:
  ret
}