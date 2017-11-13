#' Profile Parsimony Score
#'
#' Calculate a tree's Profile Parsimony score with a given dataset, after Faith and Trueman (2001)
#'
#' @template treeParam
#' @param dataset Dataset of class \code{phyDat}.  The dataset should have been prepared using
#'                \code{dataset <- \link{PrepareDataProfile}(dataset)}; if this step
#'                has not been completed, the dataset will 
#'                be (time-consumingly) prepared within the function.
#'                In subsidiary functions, the dataset will have been initialized using 
#'                \code{ProfileInitMorphy}, must be destroyed using \code{ProfileDestroyMorphy}.
#'
#' @return Zero minus the profile score (because the optimization algorithm assumes that
#'         smaller numbers are better)
#'
#' @references
#'    Faith, D. P. & Trueman, J. W. H. (2001). \cite{Towards an inclusive philosophy for 
#'    phylogenetic inference.} Systematic Biology 50:3, 331-350, doi: 
#'    \href{http://dx.doi.org/10.1080/10635150118627}{10.1080/10635150118627}
#'
#' @examples
#'   data(referenceTree)
#'   data(congreveLamsdellMatrices)
#'   # In actual use, the dataset should be prepared with a much higher precision: try 1e+06?
#'   # Of course, gaining higher precision takes substantially longer.
#'   dataset <- PrepareDataProfile(congreveLamsdellMatrices[[42]], precision=1e+03)
#'   ProfileScore(referenceTree, dataset)
#'
#' @author Martin R. Smith
#' @keywords tree
#' @export
ProfileScore <- function (tree, dataset) {
  if (class(dataset) != 'phyDat') {
    stop('Invalid dataset type; prepare dataset with PhyDat() and PrepareDataProfile().')
  }
  if (!('info.amounts' %in% names(attributes(dataset)))) dataset <- PrepareDataProfile(dataset)
  at <- attributes(dataset)
  nChar  <- at$nr # strictly, transformation series patterns; these'll be upweighted later
  weight <- at$weight
  steps <- FitchSteps(tree, dataset)
  info <- at$info.amounts
  nRowInfo <- nrow(info)
  infoPerChar <- vapply(seq_len(nChar), function (i) {
      stepRow <- max(0L, steps[i] - 1L) + 1L
      return(if (stepRow > nRowInfo) 0 else info[stepRow, i])
    }, double(1))
  # Return:
  -sum(infoPerChar * weight)
}

#' @describeIn ProfileScore Scorer for initialized dataset.
#' @template treeParent
#' @template treeChild
#' @export
ProfileScoreMorphy <- function (parent, child, dataset) {
  steps <- vapply(attr(dataset, 'morphyObjs'), MorphyLength, parent=parent, child=child, integer(1))
  info <- attr(dataset, 'info.amounts')
  nRowInfo <- nrow(info)
  # Return:
  -sum(vapply(seq_along(steps), function (i) {
    stepRow <- max(0L, steps[i] - 1L) + 1L
    return(if (stepRow > nRowInfo) 0 else info[stepRow, i])
  }, double(1)) * attr(dataset, 'weight'))
}

#' @describeIn ProfileScore Initialize dataset by adding morphyObjs.
#' @export
ProfileInitMorphy <- function (dataset) {
  attr(dataset, 'morphyObjs') <- apply(dataset, 1, SingleCharMorphy)
  # Return:
  dataset
}

#' @describeIn ProfileScore Free memory from morphyObjs initialized by \kbd{ProfileScoreMorphy}.
#' @export
ProfileDestroyMorphy <- function (dataset) {
  vapply(attr(dataset, 'morphyObjs'), UnloadMorphy, integer(1))
}

ProfileTreeSearch <- function (tree, dataset, Rearrange = RootedTBR,
                        maxIter = 100, maxHits = 20, forestSize = 1,
                        verbosity = 1, precision=40000, ...) {
  if (class(dataset) != 'phyDat') stop("Unrecognized dataset class; should be phyDat, not ", class(dataset), '.')
  if (!('info.amounts' %in% names(attributes(dataset)))) dataset <- PrepareDataProfile(dataset, precision)
  at <- attributes(dataset)
  
  TreeSearch(tree, dataset, nChar=at$nr, weight=at$weight, info=at$info.amounts,
             nRowInfo=nrow(at$info.amounts), 
             InitializeData = ProfileInitMorphy,
             CleanUpData = ProfileDestroyMorphy,
             TreeScorer = ProfileScoreMorphy,
             Rearrange = Rearrange, 
             maxIter = maxIter, maxHits = maxhits, forestSize = forestSize,
             verbosity = verbosity, ...)
}