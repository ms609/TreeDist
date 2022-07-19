.Bercow <- function (x) order(order(x, method = "radix"), method = "radix")

.Triangle <- function (n) n * (n + 1) / 2

.TrustSum <- function(r, rHat, N, k) {
  sum(vapply(seq_len(N), function(i) {
    Uk <- rHat[, i] <= k & r[, i] > k
    sum(r[Uk, i] - k, na.rm = TRUE)
  }, double(1)))
}

.MMax <- function (N, k) {
  UkMax <- min(k, N - k)
  if (UkMax < 1) {
    warning("All points are nearest neighbours. Decrease `neighbours`.")
  }
  # Return:
  N * ((UkMax * (N - k - 1)) - .Triangle(UkMax - 1))
}

#' Faithfulness of mapped distances
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
#' @param neighbours Integer specifying number of nearest neighbours to use in
#' calculation.  This should typically be small relative to the number of
#' points.
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
#' @template MRS
#' 
#' @references
#' \insertAllCited{}
#' @family tree space functions
#' @export
MappingQuality <- function(original, mapped, neighbours = 10L) {
  originalRank <- apply(as.matrix(original), 2, .Bercow) - 1
  mappedRank <- apply(as.matrix(mapped), 2, .Bercow) - 1
  diag(originalRank) <- diag(mappedRank) <- NA
  
  if (!identical(dim(originalRank), dim(mappedRank))) {
    stop("Original and mapped distances must have the same dimensions")
  }
  
  N <- dim(originalRank)[2]
  k <- neighbours
  MMax <- .MMax(N, k)
  
  trust <- 1 - (.TrustSum(originalRank, mappedRank, N, neighbours) / MMax)
  cont <- 1 - (.TrustSum(mappedRank, originalRank, N, neighbours) / MMax)
  txc <- trust * cont
  c('Trustworthiness' = trust,
    'Continuity' = cont,
    'TxC' = txc,
    'sqrtTxC' = sqrt(txc))
}

#' @rdname MappingQuality
#' @export
ProjectionQuality <- MappingQuality
