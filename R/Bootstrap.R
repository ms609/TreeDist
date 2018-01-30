#' Ratchet bootstrapper
#' 
#' @template edgeListParam
#' @template morphyObjParam
#' @template EdgeSwapperParam
#' @param maxIter maximum number of iterations to perform in tree search
#' @param maxHits maximum number of hits to accomplish in tree search
#' @template verbosityParam
#' @param \dots further parameters to send to \code{TreeScorer}
#'
#' @return A tree that is optimal under a random sampling of the original characters
#' @export
MorphyBootstrap <- function (edgeList, morphyObj, EdgeSwapper = NNISwap, 
                             maxIter, maxHits, verbosity=1L, ...) {
  startWeights <- MorphyWeights(morphyObj)[1, ]
  eachChar <- seq_along(startWeights)
  deindexedChars <- rep(eachChar, startWeights)
  resampling <- tabulate(sample(deindexedChars, replace=TRUE), length(startWeights))
  errors <- vapply(eachChar, function (i) 
         mpl_set_charac_weight(i, resampling[i], morphyObj), integer(1))
  if (any(errors)) stop ("Error resampling morphy object: ", mpl_translate_error(unique(errors[errors < 0L])))
  if (mpl_apply_tipdata(morphyObj) -> error) stop("Error applying tip data: ", mpl_translate_error(error))
  res <- EdgeListSearch(edgeList, morphyObj, EdgeSwapper=EdgeSwapper, maxIter=maxIter, maxHits=maxHits, verbosity=verbosity-1L, ...)
  errors <- vapply(eachChar, function (i) 
         mpl_set_charac_weight(i, startWeights[i], morphyObj), integer(1))
  if (any(errors)) stop ("Error resampling morphy object: ", mpl_translate_error(unique(errors[errors < 0L])))
  if (mpl_apply_tipdata(morphyObj) -> error) stop("Error applying tip data: ", mpl_translate_error(error))
  res[1:2]
}


#' @describeIn MorphyBootstrap Generic Bootstrapper for tree scoring functions that consider each character separately
#' @template datasetParam
#' @export
#' @keywords internal
CharacterwiseBootstrap <- function (edgeList, dataset, TreeScorer, EdgeSwapper = NNISwap, 
                                   maxIter, maxHits, verbosity=1L, ...) {
  att <- attributes(dataset)
  startWeights <- att[['weight']]
  eachChar <- seq_along(startWeights)
  deindexedChars <- rep(eachChar, startWeights)
  resampling <- tabulate(sample(deindexedChars, replace=TRUE), length(startWeights))
  sampled <- resampling != 0
  sampledData <- lapply(dataset, function (x) x[sampled])
  sampledAtt <- att
  sampledAtt[['weight']] <- resampling[sampled]
  sampledAtt[['index']] <- rep(seq_len(sum(sampled)), resampling[sampled])
  sampledAtt[['info.amounts']] <- att[['info.amounts']][, sampled]
  sampledAtt[['morphyObjs']] <- att[['morphyObjs']][sampled]
  attributes(sampledData) <- sampledAtt
  
  res <- EdgeListSearch(edgeList, sampledData, TreeScorer=TreeScorer,
                        EdgeSwapper=EdgeSwapper, maxIter=maxIter, maxHits=maxHits, verbosity=verbosity-1L, ...)
  
  res[1:2]
}

#' @describeIn MorphyBootstrap Bootstrapper for Profile Parsimony
#' @template datasetParam
#' @export
ProfileBootstrap <- function (edgeList, dataset, EdgeSwapper = NNISwap, 
                              maxIter, maxHits, verbosity=1L, ...) {
  CharacterwiseBootstrap(edgeList=edgeList, dataset=dataset, TreeScorer=ProfileScoreMorphy, EdgeSwapper=EdgeSwapper, 
                         maxIter=maxIter, maxHits=maxHits, verbosity=verbosity, ...)
}


#' @describeIn MorphyBootstrap Bootstrapper for Implied weighting
#' @template datasetParam
#' @export
IWBootstrap <- function (edgeList, dataset, EdgeSwapper = NNISwap, 
                              maxIter, maxHits, verbosity=1L, ...) {
  CharacterwiseBootstrap(edgeList=edgeList, dataset=dataset, TreeScorer=IWMorphy, EdgeSwapper=EdgeSwapper, 
                         maxIter=maxIter, maxHits=maxHits, verbosity=verbosity, ...)
}

