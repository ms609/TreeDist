#' @keywords internal
.calcContinuityFromDist <- function(distReference, distLowDim, kTM) {
  distReference <- as.matrix(distReference)
  distLowDim <- as.matrix(distLowDim)
  
  stopifnot(ncol(distReference) == ncol(distLowDim))
  stopifnot(nrow(distReference) == nrow(distLowDim))
  
  rankReference <- apply(
    as.matrix(distReference), 2, function(w) order(order(w)))
  rankLowDim <- apply(
    as.matrix(distLowDim), 2, function(w) order(order(w)))
  
  .calcContinuityFromRank(
    rankReference = rankReference,
    rankLowDim = rankLowDim, kTM = kTM
  )
}

#' @keywords internal
.calcContinuityFromRank <- function(rankReference, rankLowDim, kTM) {
  stopifnot(ncol(rankReference) == ncol(rankLowDim))
  stopifnot(nrow(rankReference) == nrow(rankLowDim))
  
  N <- ncol(rankReference)
  
  1 - 2/(N * kTM * (2 * N - 3 * kTM - 1)) *
    sum(vapply(seq_len(ncol(rankReference)), function(i) {
      sum((rankLowDim[, i] - kTM) * (rankReference[, i] <= kTM) *
            (rankLowDim[, i] > kTM))
    }, NA_real_))
}

#' @keywords internal
.calcTrustworthinessFromRank <- function(rankReference, rankLowDim, kTM) {
  stopifnot(ncol(rankReference) == ncol(rankLowDim))
  stopifnot(nrow(rankReference) == nrow(rankLowDim))
  
  N <- ncol(rankReference)
  
  1 - 2/(N * kTM * (2 * N - 3 * kTM - 1)) *
    sum(vapply(seq_len(ncol(rankLowDim)), function(i) {
      sum((rankReference[, i] - kTM) * (rankLowDim[, i] <= kTM) *
            (rankReference[, i] > kTM))
    }, NA_real_))
}

#' @keywords internal
#' @source Charlotte Soneson's '[dreval](https://github.com/csoneson/dreval/blob/master/R/trustworthiness.R)'
.calcTrustworthinessFromDist <- function(distReference, distLowDim, kTM) {
  distReference <- as.matrix(distReference)
  distLowDim <- as.matrix(distLowDim)
  
  stopifnot(ncol(distReference) == ncol(distLowDim))
  stopifnot(nrow(distReference) == nrow(distLowDim))
  
  rankReference <- apply(as.matrix(distReference), 2,
                         function(w) order(order(w)))
  rankLowDim <- apply(as.matrix(distLowDim), 2,
                      function(w) order(order(w)))
  
  .calcTrustworthinessFromRank(rankReference = rankReference,
                              rankLowDim = rankLowDim, kTM = kTM)
}

#' Quality of mapped distances
#' 
#' `MappingQuality()` calculates the trustworthiness and continuity
#' of mapped distances \insertCite{Venna2001,Kaski2003}{TreeDist}.
#' Trustworthiness measures, on a scale from 0--1,
#' the degree to which points that are nearby in a mapping are truly close
#' neighbours; continuity, the extent to which points that are truly nearby 
#' retain their close spatial proximity in a mapping.
#' 
#' 
#' @param original,mapped Square matrix or `dist` object containing 
#' original / mapped pairwise distances.
#' @param neighbours Number of nearest neighbours to use in calculation.
#' 
#' @return `MappingQuality()` returns a named vector of length four, 
#' containing the entries: `Trustworthiness`, `Continuity`, `TxC` 
#' (the product of these values), and `sqrtTxC` (its square root).
#' 
#' @examples
#' library('TreeTools', quietly = TRUE, warn.conflict = FALSE)
#' trees <- as.phylo(0:10, nTip = 8)
#' distances <- ClusteringInfoDistance(trees)
#' mapping <- cmdscale(distances)
#' MappingQuality(distances, dist(mapping), 4)
#' @author Wrapper for functions from Charlotte Soneson's \pkg{dreval},
#' https://github.com/csoneson/dreval/blob/master/R/trustworthiness.R
#' 
#' @references
#' \insertAllCited{}
#' @family tree space functions
#' @export
MappingQuality <- function (original, mapped, neighbours = 10L) {
  trust <- .calcTrustworthinessFromDist(original, mapped, neighbours)
  cont <- .calcContinuityFromDist(original, mapped, neighbours)
  c('Trustworthiness' = trust,
    'Continuity' = cont,
    'TxC' = trust * cont,
    'sqrtTxC' = sqrt(trust * cont))
}

#' @rdname MappingQuality
#' @export
ProjectionQuality <- MappingQuality