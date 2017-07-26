#' @useDynLib TreeSearch FITCH
#' @keywords internal
#' @export
C_Fitch_Score <- function (characters, nChar, parent, child, nEdge, weight, maxNode, nTip) {
  .Call("FITCH", as.integer(characters), as.integer(nChar),
        as.integer(parent), as.integer(child), as.integer(nEdge),
        as.double(weight), as.integer(maxNode), as.integer(nTip))[[1]]
}
#' @useDynLib TreeSearch FITCH
#' @keywords internal
#' @export
C_Fitch_Steps <- function (characters, nChar, parent, child, nEdge, weight, maxNode, nTip) {
  .Call("FITCH", as.integer(characters), as.integer(nChar),
        as.integer(parent), as.integer(child), as.integer(nEdge),
        as.double(weight), as.integer(maxNode), as.integer(nTip))[[2]]
}
#' @useDynLib TreeSearch FITCH
#' @keywords internal
#' @export
C_Fitch <- function (characters, nChar, parent, child, nEdge, weight, maxNode, nTip) {
  .Call("FITCH", as.integer(characters), as.integer(nChar),
        as.integer(parent), as.integer(child), as.integer(nEdge),
        as.double(weight), as.integer(maxNode), as.integer(nTip))
}

#' @keywords internal
#' @export
TipsAreNames <- function(data, tips) as.integer(unlist(data[tips]))
#' @keywords internal
#' @export
TipsAreRows <- function(data, tips) as.integer(data[tips, ])
#' @keywords internal
#' @export
TipsAreColumns <- function(data, tips) as.integer(data[, tips])

#' Fitch score
#' 
#' @template treeParam
#' @param data A list or matrix whose list entries, rows or columns (as specified in TipData)
#'             correspond to the character data associated with the tips of the tree. 
#'             Likely of class \code{phyDat}.
#' @param TipData function to be used to match tree tips to dataset; probably one of 
#'                \code{TipsAreNames}, \code{TipsAreRows}, or \code{TipsAreColumns}
#' @param at Attributes of the dataset (looked up automatically if not supplied)
#' @param FitchFunction function to be used to calculte parsimony score.
#' @return A vector listing the number of 'parsimony steps' calculated for each character
#' @export
Fitch <- function (tree, data, TipData = TipsAreNames, at = attributes(data),
                        FitchFunction = C_Fitch_Score) { 
  if (is.null(at$order) || at$order == "cladewise") tree <- Postorder(tree)
  treeEdge <- tree$edge
  parent <- treeEdge[, 1]
  child <- treeEdge[, 2]
  tipLabel <- tree$tip.label
  
  return(FitchFunction(
      characters = TipData(data, tipLabel), 
      nChar = at$nr,
      parent, child, 
      nEdge = length(parent), 
      weight = at$weight, 
      maxNode = parent[1], #max(parent),
      nTip = length(tipLabel)
    ) 
    
    FitchFunction(TipData(data, tipLabel), at$nr, parent, child, length(parent), at$weight, parent[1], nTip = length(tipLabel))
  )
}
#' @describeIn Fitch returns the parsimony score only
#' @export
FitchScore <- function (tree, data, TipData = TipsAreNames, at = attributes(data))
  Fitch(tree, data, TipData, at, C_Fitch_Score)

#' @describeIn Fitch returns a vector listing the number of steps for each character
#' @export
FitchSteps <- function (tree, data, TipData = TipsAreNames, at = attributes(data))
  Fitch(tree, data, TipData, at, C_Fitch_Steps)
