#' Bootstrap tree search with inapplicable data
#' 
#' @template labelledTreeParam

#'
#' @return A tree that is optimal under a random sampling of the original characters
#' @export  
BootstrapTree <- function (phy, x, maxIter, maxHits, TreeScorer = FitchScore, track=1, ...) {
## Simplified version of phangorn::bootstrap.phyDat, with bs=1 and multicore=FALSE
  at <- attributes(x)
  bootstrappedWeight <- BootstrapWeightings(at)
  keep <- bootstrappedWeight > 0
  if (is.null(at$names)) x <- x[keep, ] else x <- lapply(x, function (i) i[keep])
  xDim <- dim(x)
  xDimNames <- dimnames(x)
  mostattributes(x) <- at
  dim(x) <- xDim
  dimnames(x) <- xDimNames
  attr(x, 'weight') <- bootstrappedWeight[keep]
  for (attrName in at$bootstrap) {
    attribute <- at[attrName][[1]]
    if (length(dim(attribute)) == 2) {
      if (dim(attribute)[2] != length(keep)) warning("Bootstrap mismatch for attribute", attrName)
      attr(x, attrName) <- attribute[, keep]
    } else if (length(attribute) == length(keep)) {
      attr(x, attrName) <- attribute[keep]
    } 
  }
  attr(x, 'nr') <- sum(keep)
  
  attr(phy, 'score') <- NULL
  res <- DoTreeSearch(phy, x, TreeScorer=TreeScorer, method='NNI', maxIter=maxIter,
                      maxHits=maxHits, track=track-1, ...)
  attr(res, 'score') <- NULL
  attr(res, 'hits') <- NULL
  res
}

#' @keywords internal
#' @export
BootstrapWeightings <- function (attribs) {
  weight <- attribs$weight
  v <- rep(1:length(weight), weight)
  BS <- tabulate(sample(v, replace=TRUE), length(weight)) 
}