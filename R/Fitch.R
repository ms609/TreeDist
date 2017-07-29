#' Fitch Score
#' Calculate the parsimony score of a tree (number of steps) using the Fitch algoritgh
#' @return the parsimony score (an integer)
#' @useDynLib TreeSearch FITCH
#' @keywords internal
#' @export
C_Fitch_Score <- function (characters, nChar, parent, child, nEdge, weight, maxNode, nTip) {
  .Call("FITCH", as.integer(characters), as.integer(nChar),
        as.integer(parent), as.integer(child), as.integer(nEdge),
        as.double(weight), as.integer(maxNode), as.integer(nTip))[[1]]
}
#' Fitch steps
#' @return the number of steps the tree enforces on each character
#' @useDynLib TreeSearch FITCH
#' @keywords internal
#' @export
C_Fitch_Steps <- function (characters, nChar, parent, child, nEdge, weight, maxNode, nTip) {
  .Call("FITCH", as.integer(characters), as.integer(nChar),
        as.integer(parent), as.integer(child), as.integer(nEdge),
        as.double(weight), as.integer(maxNode), as.integer(nTip))[[2]]
}
#' Wrapper to FITCH
#' @return the full return of the phangorn C function FITCH
#' @useDynLib TreeSearch FITCH
#' @keywords internal
#' @export
C_Fitch <- function (characters, nChar, parent, child, nEdge, weight, maxNode, nTip) {
  .Call("FITCH", as.integer(characters), as.integer(nChar),
        as.integer(parent), as.integer(child), as.integer(nEdge),
        as.double(weight), as.integer(maxNode), as.integer(nTip))
}

#' Extract character data from dataset
#'
#' Speficies how characters are stored in a data object
#' 
#' @param data a matrix or list containing symbols associated with each tip; for example,
#'             a dataset of class \code{phyDat}
#' @param tips vector detailing the tips to be selected, whether as their
#'        names or numbers corresponding to their rows/columns
#' @return an integer vector, listing the tokens associated with each character for each tip in turn
#'         - ready to send to FITCH or equivalent C routine
#' @keywords internal
#' @export
TipsAreNames <- function(data, tips) as.integer(unlist(data[tips]))
#' @describeIn TipsAreNames use if each row in a matrix corresponds to a tip
#' @keywords internal
#' @export
TipsAreRows <- function(data, tips) as.integer(data[tips, ])
#' @describeIn TipsAreNames use if each column in a matrix corresponds to a tip
#' @keywords internal
#' @export
TipsAreColumns <- function(data, tips) as.integer(data[, tips])

#' Fitch score
#' 
#' @template treeParam
#' @param data A list or matrix whose list entries, rows or columns (as specified in TipData)
#'             correspond to the character data associated with the tips of the tree. 
#'             Likely of class \code{phyDat}.  The number of characters should be stored
#'             as the attribute data$nr, and the weight of each character in data$weight;
#'             if data is a matrix, they will be calculated automatically (weights set to 1).
#' @param TipData function to be used to match tree tips to dataset; probably one of 
#'                \code{TipsAreNames}, \code{TipsAreRows}, or \code{TipsAreColumns}
#' @param at Attributes of the dataset (looked up automatically if not supplied)
#' @param FitchFunction function to be used to calculte parsimony score.
#' @return A vector listing the number of 'parsimony steps' calculated for each character
#' @export
Fitch <- function (tree, data, TipData = TipsAreNames, at = attributes(data),
                        FitchFunction = C_Fitch_Score) { 
  treeOrder <- attr(tree, 'order')
  if (is.null(treeOrder) || treeOrder == "cladewise") tree <- Postorder(tree)
  treeEdge <- tree$edge
  parent <- treeEdge[, 1]
  child <- treeEdge[, 2]
  tipLabel <- tree$tip.label
  nr <- at$nr
  if (is.null(at$nr)) {
    nr <- dim(data)
    nr <- if (nr[1] == length(tipLabel)) nr[2] else nr[1]
  }
  charWeights <- at$weight
  if (is.null(charWeights)) charWeights <- rep(1, nr)
  if (class(data) =='phyDat') {
    levs <- attr(data, 'levels')
    contrast <- attr(data, 'contrast')
    index <- as.integer(contrast %*% 2L ^ (seq_along(attr(data, 'levels')) - 1))
    data <- vapply(data, function (X) index[X], integer(nr))
    characters <- TipsAreColumns(data, tipLabel)
  } else {
    characters <- TipData(data, tipLabel)
  }
  
  return(FitchFunction(
      characters = characters, 
      nChar = nr,
      parent, child,
      nEdge = length(parent), 
      weight = charWeights, 
      maxNode = parent[1], #max(parent),
      nTip = length(tipLabel)
    )
  )
  # TODO DELTE debugging line:
  #    FitchFunction(TipData(data, tipLabel), at$nr, parent, child, length(parent), at$weight, parent[1], nTip = length(tipLabel))
}
#' @describeIn Fitch returns the parsimony score only
#' @export
FitchScore <- function (tree, data, TipData = NULL, at = attributes(data)) {
  if (is.null(TipData)) {
    if (class(data) != 'phyDat') stop("Expected a phyDat object. Please specify the TipData function relevant
      to the class of data that you have provided")
    TipData <- TipsAreNames
    if (!all(tree$tip.label %in% names(data))) stop ("Some tips could not be found in the data provided.")
  }
  Fitch(tree, data, TipData, at, C_Fitch_Score)
}

#' @describeIn Fitch returns the parsimony score only, without checking that data is well formatted
#' @keywords internal
#' @export
FasterFitchScore <- function (tree, data, TipData = TipsAreNames, at = attributes(data))
  Fitch(tree, data, TipData, at, C_Fitch_Score)

  
#' @describeIn Fitch returns a vector listing the number of steps for each character
#' @export
FitchSteps <- function (tree, data, TipData = TipsAreNames, at = attributes(data))
  Fitch(tree, data, TipData, at, C_Fitch_Steps)
