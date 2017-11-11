#' Fitch Score
#' Calculate the parsimony score of a tree (number of steps) using the Fitch algorithm
#' @return the parsimony score (an integer)
#' @keywords internal
#' @export
C_Fitch_Score <- function (characters, nChar, parent, child, nEdge, weight, maxNode, nTip) {
  C_Fitch(characters, nChar, parent, child, nEdge, weight, maxNode, nTip)[[1]]
}
#' Fitch steps
#' @return the number of steps the tree enforces on each character
#' @keywords internal
#' @export
C_Fitch_Steps <- function (characters, nChar, parent, child, nEdge, weight, maxNode, nTip) {
  C_Fitch(characters, nChar, parent, child, nEdge, weight, maxNode, nTip)[[2]]
}
#' Wrapper to FITCH
#' @return the full return of the phangorn C function FITCH
#' @useDynLib TreeSearch, .registration = TRUE, .fixes="C_"
#' @keywords internal
#' @export
C_Fitch <- function (characters, nChar, parent, child, nEdge, weight, maxNode, nTip) {
  .Call(C_phangorn_FITCH, as.integer(characters), as.integer(nChar),
        as.integer(parent), as.integer(child), as.integer(nEdge),
        as.double(weight), as.integer(maxNode), as.integer(nTip))
}
###   
###   #' @describeIn C_Fitch Checks parameters before running \code{C_Fitch}
###   #' @keywords internal
###   #' @export
###   C_Fitch_Checks <- function (characters, nChar, parent, child, nEdge, weight, maxNode, nTip) {
###     iCharacters <- as.integer(characters)
###     Assert(length(iCharacters) == nChar * nTip)
###     Assert(length(parent) == length(child))
###     Assert(length(parent) == nEdge)
###     Assert(length(weight) == nChar)
###     Assert(maxNode == max(parent))
###     Assert(parent[1] == max(parent)) # Are nodes numbered in Postorder??    
###     C_Fitch(characters, nChar, parent, child, nEdge, weight, maxNode, nTip)
###   }

#' Extract character data from dataset
#'
#' Specifies how characters are stored in a dataset object
#' 
#' @param dataset a matrix or list containing symbols associated with each tip; for example,
#'             a dataset of class \code{phyDat}
#' @param tips vector detailing the tips to be selected, whether as their
#'        names or numbers corresponding to their rows/columns
#' @return an integer vector, listing the tokens associated with each character for each tip in turn
#'         - ready to send to FITCH or equivalent C routine
#' @keywords internal
#' @export
TipsAreNames <- function(dataset, tips) as.integer(unlist(dataset[tips]))

###   #' @describeIn TipsAreNames use if each row in a matrix corresponds to a tip
###   #' @keywords internal
###   #' @export
###   TipsAreRows <- function(dataset, tips) as.integer(dataset[tips, ])

#' @describeIn TipsAreNames use if each column in a matrix corresponds to a tip
#' @keywords internal
#' @export
TipsAreColumns <- function(dataset, tips) as.integer(dataset[, tips])

#' Fitch score
#' 
#' @template treeParam
#' @param dataset A list or matrix whose list entries, rows or columns (as specified in TipData)
#'             correspond to the character data associated with the tips of the tree. 
#'             Likely of class \code{phyDat}.  The number of characters should be stored
#'             as the attribute dataset$nr, and the weight of each character in dataset$weight;
#'             if dataset is a matrix, they will be calculated automatically (weights set to 1).
#' @param TipData function to be used to match tree tips to dataset; probably one of 
#'                \code{TipsAreNames}, \code{TipsAreRows}, or \code{TipsAreColumns}
#' @param at Attributes of the dataset (looked up automatically if not supplied)
#' @param FitchFunction function to be used to calculate parsimony score.
#'
#' @return A vector listing the number of 'parsimony steps' calculated for each character
#'
#' @importFrom phangorn fitch
#' @export
Fitch <- function (tree, dataset, TipData = TipsAreNames, at = attributes(dataset),
                        FitchFunction = C_Fitch_Score) {
  treeOrder <- attr(tree, 'order')
  if (is.null(treeOrder) || treeOrder != 'postorder') tree <- Postorder(tree)
  treeEdge <- tree$edge
  parent <- treeEdge[, 1]
  child <- treeEdge[, 2]
  tipLabel <- tree$tip.label
  nr <- at$nr
  if (is.null(at$nr)) {
    nr <- dim(dataset)
    nr <- if (nr[1] == length(tipLabel)) nr[2] else nr[1]
  }
  charWeights <- at$weight
  if (is.null(charWeights)) charWeights <- rep(1, nr)
  if (class(dataset) =='phyDat') {
    levs <- attr(dataset, 'levels')
    contrast <- attr(dataset, 'contrast')
    index <- as.integer(contrast %*% 2L ^ (seq_along(attr(dataset, 'levels')) - 1))
    reformedData <- vapply(dataset, function (X) index[X], integer(nr))
    characters <- TipsAreColumns(reformedData, tipLabel)
  } else {
    characters <- TipData(dataset, tipLabel)
  }
  
  return(FitchFunction(
      characters = characters, 
      nChar = nr,
      parent, child,
      nEdge = length(parent), 
      weight = charWeights, 
      maxNode = max(parent), #parent[1] IF tree in postorder
      nTip = length(tipLabel)
    )
  )
  # TODO DELTE debugging line:
  #    FitchFunction(TipData(dataset, tipLabel), at$nr, parent, child, length(parent), at$weight, parent[1], nTip = length(tipLabel))
}

#' @describeIn Fitch returns the parsimony score only
#' @export
FitchScore <- function (tree, dataset, TipData = NULL, at = attributes(dataset)) {
  if (is.null(TipData)) {
    if (class(dataset) != 'phyDat') stop("Expected a phyDat object. Please specify the TipData function relevant
      to the class of dataset that you have provided")
    TipData <- TipsAreNames
    if (!all(tree$tip.label %in% names(dataset))) stop ("Some tips could not be found in the dataset provided.")
  }
  Fitch(tree, dataset, TipData, at, C_Fitch_Score)
}

###   #' @describeIn Fitch returns the parsimony score only, when a Fitch instance has already been initiated
###   #' @export
###   IFitchScore <- function (tree, nChar) {
###     phangorn:::fast.fitch(tree, nChar, ps=FALSE)
###   }
###   
###   #' @describeIn Fitch returns the parsimony score only, without checking that dataset is well formatted
###   #' @keywords internal
###   #' @export
###   FasterFitchScore <- function (tree, dataset, TipData = TipsAreNames, at = attributes(dataset))
###     Fitch(tree, dataset, TipData, at, C_Fitch_Score)
###
###
###   IFitchSteps <- function (tree, nChar) {
###     phangorn:::fast.fitch(tree, nChar, ps=FALSE)
###   }
###     
###   #' @describeIn Fitch returns a vector listing the number of steps for each character
###   #' @export
###   FitchSteps <- function (tree, dataset, TipData = TipsAreNames, at = attributes(dataset))
###     Fitch(tree, dataset, TipData, at, C_Fitch_Steps)
