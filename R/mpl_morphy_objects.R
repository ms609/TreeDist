#' Details the attributes of a morphy object
#'
#' @param object A Morphy object
#' @param \dots any other parameters... 
#'
#' @return A list detailing the number of taxa, internal nodes, and characters and their weights.
#'
#' @author Martin R. Smith
#' @method summary morphyPtr
#' @export
summary.morphyPtr <- function (object, ...) {
  ans <- list()
  class(ans) <- "summary.morphyPtr"
  nTax <- mpl_get_numtaxa(object)
  nChar <- mpl_get_num_charac(object)
  charWeights <- MorphyWeights(object)

  ans$nTax <- nTax 
  ans$nChar <- nChar 
  ans$nInternal <- mpl_get_num_internal_nodes(object)
  ans$charWeights <- charWeights
  ans$allStates <- mpl_get_symbols(object)
  # Return:
  ans
}

#' Report the character weightings associated with a Morphy object
#'
#' @template morphyObjParam
#' @return a matrix of dimensions (2, number of characters); row 1 lists the
#'         exact rates specified by the user; row 2 the approximate (integral)
#'         weights used by MorphyLib
#'
#' @author Martin R. Smith
#' @export
MorphyWeights <- function(morphyObj) {
 charWeights <- vapply(seq_len(mpl_get_num_charac(morphyObj)), mpl_get_charac_weight, list(0, 0), morphyobj=morphyObj)
 dimnames(charWeights) <- list(c("exact", "approx"), NULL)
 charWeights
}

#' Set the character weightings associated with a Morphy object
#'
#' @param weight A vector listing the new weights to be applied to each character
#' @template morphyObjParam
#' @param checkInput Whether to sanity-check input data before applying. 
#'         Defaults to TRUE to protect the user from crashes.
#'
#' @return The Morphy error code generated when applying tipData
#' 
#' @author Martin R. Smith
#' @export
SetMorphyWeights <- function (weight, morphyObj, checkInput = TRUE) {
  if (checkInput) if (length(weight) != mpl_get_num_charac(morphyObj)) stop("Number of weights not equal to number of character entries")
  errors <- vapply(seq_along(weight), function (i)
    mpl_set_charac_weight(i, weight[i], morphyObj), integer(1))
  if(any(errors != 0)) warning("Morphy Error encountered: ", mpl_translate_error(errors[errors<0]))
  mpl_apply_tipdata(morphyObj)
}

#' Initialize a Morphy Object from a phyDat object
#' 
#' Creates a new Morphy object with the same size and characters as the phyDat object 
#'
#' @param phy An object of class \code{\link{phyDat}}.
#' @return A pointer to an initialized Morphy object.
#' 
#' @author Martin R. Smith
#' @importFrom phangorn phyDat
#' @export
PhyDat2Morphy <- function (phy) {
  
  if (class(phy) != 'phyDat') {
    stop('Invalid data type ', class(phy), '; should be phyDat.')
  }
  
  morphyObj <- mpl_new_Morphy()
  nTax <- length(phy)
  weight <- attr(phy, 'weight')
  nChar <- attr(phy, 'nr')
  
  if (mpl_init_Morphy(nTax, nChar, morphyObj) -> error) {
    stop("Error ", mpl_translate_error(error), " in mpl_init_Morphy")
  }
  if (mpl_attach_rawdata(PhyToString(phy, ps=';', useIndex=FALSE, byTaxon=TRUE, concatenate=TRUE),
                          morphyObj) -> error) {
    stop("Error ", mpl_translate_error(error), " in mpl_attach_rawdata")
  }
  if (mpl_set_num_internal_nodes(nTax - 1L, morphyObj) -> error) { # One is the 'dummy root'
    stop("Error ", mpl_translate_error(error), " in mpl_set_num_internal_nodes")
  }
  if (any(vapply(seq_len(nChar), function (i) mpl_set_parsim_t(i, 'FITCH', morphyObj), integer(1)) 
      -> error)) {
      stop("Error ", mpl_translate_error(min(error)), "in mpl_set_parsim_t")
  }
  if (any(vapply(seq_len(nChar), function (i) mpl_set_charac_weight(i, weight[i], morphyObj),
      integer(1)) -> error)) {
    stop("Error ", mpl_translate_error(min(error)), "in mpl_set_charac_weight")
  }
  if (mpl_apply_tipdata(morphyObj) -> error) {
    stop("Error ", mpl_translate_error(error), "in mpl_apply_tipdata")
  }
  class(morphyObj) <- 'morphyPtr'
  # Return:
  morphyObj
}

#' Check for error whilst modifying Morphy object
#' @param action action to perform
#'
#' @keywords internal
#' @export
MorphyErrorCheck <- function (action) {
  if (action -> error) stop("Morphy object encountered error ", mpl_translate_error(error), "\n")
}

#' Morphy object from single character
#' 
#' @param char State of each character at each tip in turn, in a format that will be converted
#'             to a character string by \code{\link{paste0}(char, ';', collapse='')}.
#'
#' @return A pointer to a morphyObj. Don't forget to unload it when you're done with it:
#'         \code{morphyObj <- \link{UnloadMorphy}(morphyObj)}.
#'
#' @author Martin R. Smith
#'
#' @export
SingleCharMorphy <- function (char) {
  char <- paste0(c(char, ';'), collapse='')
  entries <- gregexpr("\\{[^\\{]+\\}|\\([^\\()]+\\)|[^;]", char)
  nTip <- length(entries[[1]])
  morphyObj <- mpl_new_Morphy()
  MorphyErrorCheck(mpl_init_Morphy(nTip, 1, morphyObj)) 
  MorphyErrorCheck(mpl_attach_rawdata(char, morphyObj)) 
  MorphyErrorCheck(mpl_set_num_internal_nodes(nTip - 1L, morphyObj)) 
  MorphyErrorCheck(mpl_set_parsim_t(1, 'FITCH', morphyObj))
  MorphyErrorCheck(mpl_set_charac_weight(1, 1, morphyObj)) 
  MorphyErrorCheck(mpl_apply_tipdata(morphyObj))
  class(morphyObj) <- 'morphyPtr'
  morphyObj
}

#' Destroy a Morphy Object
#'
#' Best practice is to call \code{morphyObj <- UnloadMorphy(morphyObj)}
#' Failure to do so will cause a crash if UnloadMorphy is called on an object that 
#' has already been destroyed
#'
#' @template morphyObjParam
#' @return Morphy error code, decipherable using \code{\link{mpl_translate_error}}
#' @author Martin R. Smith
#' @export
UnloadMorphy <- function (morphyObj) {
  if (class(morphyObj) != 'morphyPtr') stop ("Object is not a valid pointer; cannot destroy.")
  if (mpl_delete_Morphy(morphyObj) -> error) {
    stop("Error ", mpl_translate_error(error), "in mpl_delete_Morphy")
  }
  return (error)
}

#' Unload this library
#' @export
UnloadTreeSearch <- function () detach("package:TreeSearch", unload=TRUE)
