#' Bootstrap tree search with inapplicable data
#' 
#' @template preorderTreeParam
#' @template datasetTreeScorerParams
#' @param maxIter maximum number of iterations to perform in tree search
#' @param maxHits number of times to find optimal tree length before stopping tree search
#' @param rooted set to FALSE if position of root may change, TRUE if position and composition of
#'               outgroup is fixed
#' @template verbosityParam
#' @template treeScorerDots
#'
#' @return A tree that is optimal under a random sampling of the original characters
#' @export  
BootstrapTree <- function (tree, dataset, TreeScorer = FitchScore, 
                           maxIter, maxHits, rooted = TRUE, verbosity=1, ...) {
## Simplified version of phangorn::bootstrap.phyDat, with bs=1 and multicore=FALSE
  at <- attributes(dataset)
  bootstrappedWeight <- BootstrapWeightings(at)
  keep <- bootstrappedWeight > 0
  if (is.null(at$names)) dataset <- dataset[keep, ] else dataset <- lapply(dataset, function (i) i[keep])
  xDim <- dim(dataset)
  xDimNames <- dimnames(dataset)
  mostattributes(dataset) <- at
  dim(dataset) <- xDim
  dimnames(dataset) <- xDimNames
  names(dataset) <- at$names
  attr(dataset, 'weight') <- bootstrappedWeight[keep]
  for (attrName in at$bootstrap) {
    attribute <- at[attrName][[1]]
    if (length(dim(attribute)) == 2) {
      if (dim(attribute)[2] != length(keep)) warning("Bootstrap mismatch for attribute", attrName)
      attr(dataset, attrName) <- attribute[, keep]
    } else if (length(attribute) == length(keep)) {
      attr(dataset, attrName) <- attribute[keep]
    }
  }
  attr(dataset, 'nr') <- sum(keep)
  ## Not sure that this is a good idea...
  ## atNames <- names(at)
  ## for (attrName in atNames[!atNames %in% names(attributes(dataset))]) {
  ##   attr(dataset, attrName) <- at[[attrName]]
  ## }
  
  if (is.null(treeOrder <- attr(tree, 'order')) || treeOrder != 'preorder') tree <- Preorder(tree)
  attr(tree, 'score') <- NULL
  res <- TreeSearch(tree, dataset, TreeScorer=TreeScorer, Rearrange=if (rooted) RootedNNI else NNI,
                      maxIter=maxIter, maxHits=maxHits, verbosity=max(0, verbosity - 1), ...)
  attr(res, 'score') <- NULL
  attr(res, 'hits') <- NULL
  res
}

#' Bootstrap Weightings
#' 
#' Resamples edges according to their prior weights
#' @param attribs Attribtues of a phyDat object
#' @return a vector giving bootstrapped weights for each character, corresponding to attribs$weights
#' @keywords internal
#' @export
BootstrapWeightings <- function (attribs) {
  weight <- attribs$weight
  v <- rep(1:length(weight), weight)
  BS <- tabulate(sample(v, replace=TRUE), length(weight)) 
}