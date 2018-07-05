#' @describeIn Ratchet Jackknife resampling. Note that at present this assumes that 
#' `InitializeData` will return a morphy object; if this doesn't hold for you, please
#' let me know and I'll make the function more general.
#' @template EdgeSwapperParam
#' @param resampleFreq Double between 0 and 1 stating proportion of characters to resample
#' @param jackIter Integer specifying number of jackknife iterations to conduct
#' @return a list of trees recovered after jackknife iterations
#' @author Martin R. Smith
#' @export
Jackknife <- function (tree, dataset, resampleFreq = 2/3,
                       InitializeData = PhyDat2Morphy,
                       CleanUpData    = UnloadMorphy,
                       TreeScorer     = MorphyLength,
                       EdgeSwapper    = TBRSwap,
                       jackIter = 5000,
                       searchIter = 4000, searchHits = 42,
                       verbosity = 1L, ...) {
  # initialize tree and data
  if (dim(tree$edge)[1] != 2 * tree$Nnode) stop("tree must be bifurcating; try rooting with ape::root")
  tree <- RenumberTips(tree, names(dataset))
  edgeList <- MatrixToList(tree$edge)
  edgeList <- RenumberEdges(edgeList[[1]], edgeList[[2]])
  
  morphyObj <- InitializeData(dataset)
  on.exit(morphyObj <- CleanUpData(morphyObj))
  
  startWeights <- MorphyWeights(morphyObj)[1, ]
  eachChar <- seq_along(startWeights)
  deindexedChars <- rep(eachChar, startWeights)
  charsToKeep <- ceiling(resampleFreq * length(deindexedChars))
  if (charsToKeep < 1L) {
    stop("resampleFreq of ", resampleFreq, " is too low; can't keep 0 of ",
         length(deindexedChars), " characters.")
  } else if (charsToKeep >= length(deindexedChars)) {
    stop("resampleFreq of ", resampleFreq, " is too high; can't keep all ",
         length(deindexedChars), " characters.")
  }
  if (verbosity > 10L) cat("\n * Beginning search:")
  
  # Conduct jackIter replicates:
  jackEdges <- vapply(seq_len(jackIter), function (x) {
    if (verbosity > 0L) {
      cat("\n * Jackknife iteration", x, "/", jackIter)
    }
    resampling <- tabulate(sample(deindexedChars, charsToKeep, replace=FALSE),
                           nbins=length(startWeights))
    errors <- vapply(eachChar, function (i) 
      mpl_set_charac_weight(i, resampling[i], morphyObj), integer(1))
    if (any(errors)) {
      stop ("Error resampling morphy object: ", mpl_translate_error(unique(errors[errors < 0L])))
    }
    if (mpl_apply_tipdata(morphyObj) -> error) {
      stop("Error applying tip data: ", mpl_translate_error(error))
    }
    res <- EdgeListSearch(edgeList[1:2], morphyObj, EdgeSwapper=EdgeSwapper,
                          maxIter=searchIter, maxHits=searchHits,
                          verbosity=verbosity-1L, ...)
    res[1:2]
  }, edgeList)
  
  jackTrees <- apply(jackEdges, 2, function(edgeList) {
    ret <- tree
    ret$edge <- ListToMatrix(edgeList)
    ret
  })
  class(jackTrees) <- 'multiPhylo'
  jackTrees
}
