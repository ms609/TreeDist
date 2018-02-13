#' Implied weights parsimony Score
#'
#' Calculate a tree's Parsimony score with a given dataset using implied weights (Goloboff 1997)
#'
#' @template treeParam
#' @param dataset Dataset of class \code{phyDat}.  The dataset should have been prepared using
#'                \code{dataset <- \link{PrepareDataIW}(dataset)}; if this step
#'                has not been completed, the dataset will 
#'                be (time-consumingly) prepared within the function.
#'                In subsidiary functions, the dataset will have been initialized using 
#'                \code{IWInitMorphy}, must be destroyed using \code{IWDestroyMorphy}.
#' @template concavityParam
#'
#' @return The 'fit', `h / h + k`, where `h` is the amount of homoplasy ('extra steps') 
#'         and `k` is a constant (the 'concavity constant')
#'
#' @references
#'  \insertRef{Goloboff1997}{TreeSearch}
#'
#' @examples
#'   data(referenceTree)
#'   data(congreveLamsdellMatrices)
#'   dataset <- PrepareDataIW(congreveLamsdellMatrices[[42]])
#'   IWScore(referenceTree, dataset)
#'
#' @author Martin R. Smith
#' @keywords tree
#' @export
IWScore <- function (tree, dataset, concavity=4) {
  if (class(dataset) != 'phyDat') {
    stop('Invalid dataset type; prepare dataset with PhyDat() and PrepareDataIW().')
  }
  if (!('info.amounts' %in% names(attributes(dataset)))) dataset <- PrepareDataIW(dataset)
  at <- attributes(dataset)
  nChar  <- at$nr # strictly, transformation series patterns; these'll be upweighted later
  weight <- at$weight
  steps <- FitchSteps(tree, dataset)
  minSteps <- at$min.steps
  homoplasies <- steps - minSteps
  
  # TODO remove this paranoid check once 100% happy with minSteps calculation.
  if (any(homoplasies < 0)) stop("Minimum steps have been miscalculated.\n", 
                            "       Please report this bug at:\n", 
                            "       https://github.com/ms609/TreeSearch/issues/new")
  fit <- homoplasies / (homoplasies + concavity)
  # Return:
  sum(fit * weight)
}

#' @describeIn ProfileScore Scorer for initialized dataset.
#' @template concavityParam
#' @param minSteps Integer vector specifying the minimum number of steps
#'                 possible for each character in `dataset`, perhaps calculated
#'                 using \code{\link{MinimumSteps}}.
#' @export
IWScoreMorphy <- function (parent, child, dataset, concavity=4, minSteps = attr(dataset, 'min.steps')) {
  steps <- vapply(attr(dataset, 'morphyObjs'), MorphyLength, parent=parent, child=child, integer(1))
  homoplasies <- steps - minSteps
  fit <- homoplasies / (homoplasies + concavity)
  # Return:
  sum(fit * attr(dataset, 'weight'))
}

#' @describeIn IWScore Initialize dataset by adding morphyObjs and min.steps.
#' @export
IWInitMorphy <- function (dataset) {
  attr(dataset, 'morphyObjs') <- 
    lapply(PhyToString(dataset, byTaxon=FALSE, useIndex=FALSE, concatenate=FALSE), 
           SingleCharMorphy)
  
  # Return:
  dataset
}


#' @describeIn TreeSearch Search using profile parsimony
#' @template concavityParam
#' @export
IWTreeSearch <- function (tree, dataset, concavity = 4, EdgeSwapper = RootedTBR,
                        maxIter = 100, maxHits = 20, forestSize = 1,
                        verbosity = 1, ...) {
  if (class(dataset) != 'phyDat') stop("Unrecognized dataset class; should be phyDat, not ", class(dataset), '.')
  if (!('min.steps' %in% names(attributes(dataset)))) dataset <- PrepareDataIW(dataset)
  at <- attributes(dataset)
  
  TreeSearch(tree, dataset, nChar=at$nr, weight=at$weight, 
             minSteps=at$min.steps, concavity = concavity,
             InitializeData = IWInitMorphy,
             CleanUpData = IWDestroyMorphy,
             TreeScorer = IWScoreMorphy,
             EdgeSwapper = EdgeSwapper, 
             maxIter = maxIter, maxHits = maxHits, forestSize = forestSize,
             verbosity = verbosity, ...)
}