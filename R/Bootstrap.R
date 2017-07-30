#' Bootstrap tree search with inapplicable data
#' 
#' @template preorderTreeParam
#' @param x a dataset in the format required by TreeScorer
#' @param maxIter maximum number of iterations to perform in tree search
#' @param maxHits number of times to find optimal tree length before stopping tree search
#' @template TreeScorerParam
#' @param rooted set to FALSE if position of root may change, TRUE if position and composition of
#'               outgroup is fixed
#' @template verbosityParam
#'
#' @return A tree that is optimal under a random sampling of the original characters
#' @export  
BootstrapTree <- function (tree, x, maxIter, maxHits, TreeScorer = FitchScore, 
  rooted = TRUE, verbosity=1, ...) {
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
  names(x) <- at$names
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
  ## Not sure that this is a good idea...
  ## atNames <- names(at)
  ## for (attrName in atNames[!atNames %in% names(attributes(x))]) {
  ##   attr(x, attrName) <- at[[attrName]]
  ## }
  
  if (is.null(treeOrder <- attr(tree, 'order')) || treeOrder != 'preorder') tree <- Preorder(tree)
  attr(tree, 'score') <- NULL
  res <- DoTreeSearch(tree, x, TreeScorer=TreeScorer, Rearrange=if (rooted) RootedNNI else NNI,
                      maxIter=maxIter, maxHits=maxHits, verbosity=max(0, verbosity - 1), ...)
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