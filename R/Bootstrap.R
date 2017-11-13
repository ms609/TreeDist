#' Bootstrap tree search with inapplicable data
#' 
#' @template edgeListParam
#' @template morphyObjParam
#' @template EdgeSwapperParam
#' @param maxIter maximum number of iterations to perform in tree search
#' @param maxHits maximum number of hits to accomplish in tree search
#' @template verbosityParam
#' @param \dots further parameters to send to \code{DoMorphySearch}
#'
#' @return A tree that is optimal under a random sampling of the original characters
#' @export
MorphyBootstrap <- function (edgeList, morphyObj, EdgeSwapper = NNISwap, 
                             maxIter, maxHits, verbosity=1L, ...) {
## Simplified version of phangorn::bootstrap.phyDat, with bs=1 and multicore=FALSE
  startWeights <- MorphyWeights(morphyObj)[1, ]
  eachChar <- seq_along(startWeights)
  v <- rep(eachChar, startWeights)
  BS <- tabulate(sample(v, replace=TRUE), length(startWeights))
  vapply(eachChar, function (i) 
         mpl_set_charac_weight(i, BS[i], morphyObj), integer(1))
  mpl_apply_tipdata(morphyObj)
  res <- DoMorphySearch(edgeList, morphyObj, EdgeSwapper=EdgeSwapper, maxIter=maxIter, maxHits=maxHits, verbosity=verbosity-1L, ...)
  vapply(eachChar, function (i) 
         mpl_set_charac_weight(i, startWeights[i], morphyObj), integer(1))
  mpl_apply_tipdata(morphyObj)
  res[1:2]
}