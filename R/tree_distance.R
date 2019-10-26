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
  splits1 <- t(as.logical(splits1)) # Convert to old-style logical matrix
  splits2 <- t(as.logical(splits2)) # Convert to old-style logical matrix
  dimSplits1 <- dim(splits1)
  dimSplits2 <- dim(splits2)
  nSplits1 <- dimSplits1[2]
  nSplits2 <- dimSplits2[2]
  if (nSplits1 == 0 || nSplits2 == 0) return (0L)
  nLeaves <- dimSplits1[1]
  if (dimSplits2[1] != nLeaves) {
    stop("Split rows must bear identical labels")
  }
  
  leafNames1 <- rownames(splits1)
  leafNames2 <- rownames(splits2)
  
  if (!is.null(leafNames2)) {
    splits2 <- unname(splits2[leafNames1, , drop=FALSE])
    splits1 <- unname(splits1) # split2[split1] faster without names
  }
  
  pairScores <- PairScorer(splits1, splits2, nSplits1, nSplits2, ...)
  
  # Return:
  TreeDistanceReturn(pairScores, maximize, reportMatching, 
                     splits1, splits2, leafNames1)
}

#' @describeIn GeneralizedRF C implementation #TODO describe
#' Note that no checks will be made to confirm that splits1 and splits2 contain
#' the same tips in the same order.  This is the responsibility of the calling
#' function.
#' @param nTip Integer specifying the number of tips in each split.
#' @references \insertRef{Jonker1987}{TreeDist}
CGRF <- function (splits1, splits2, nTip, PairScorer, 
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
      matching <- matching[seq_len(ncol(splits1))]
    }
    attr(ret, 'matching') <- matching
    
    # We're not worried about performance any more
    pairScores <- vapply(seq_len(nSplits2), function (i) {
      vapply(seq_len(nSplits1), function (j) {
        PairScorer(splits1[j, , drop = FALSE], splits2[i, , drop = FALSE],
                   nTip = nTip, ...)$score
      }, double(1))
      }, double(nSplits1))
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
#' @param leafNames Character vector listing the names corresponding to each
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
                                leafNames = NULL) {
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
    
    if (!is.null(leafNames)) {
      attr(ret, 'matchedSplits') <- 
        ReportMatching(splits1[[matched1]], splits2[[matched2]],
                       realMatch = if (maximize) {
                         pairScores[matrix(c(matched1, matched2), ncol=2L)] > 0
                       } else TRUE)
    }
  }
  # Return:
  ret
}
