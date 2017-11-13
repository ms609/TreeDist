#' @title String to phyDat
#'
#' @description Converts a PhyDat object to allow processing by MorphyDat
#'
#' @param string a string of tokens, optionally containing whitespace, with no terminating semi-colon.  Polytomies not (yet) supported; each character must correspond to a unique state, ?, or the inapplicable token (-)
#' @param tips, a character vector corresponding to the names (in order) of each taxon in the matrix
#' @param byTaxon = TRUE, string is one TAXON's coding at a time; FALSE: one CHARACTER's coding at a time
#' 
#' @examples
#' morphy <- StringToPhyDat("-?01231230?-", c('Lion', 'Gazelle'), byTaxon=TRUE)
#' # encodes the following matrix:
#' # Lion     -?0123
#' # Gazelle  1230?-
#' 
#' @template returnPhydat
#' @seealso \code{\link{phyDat}}
#' 
#' @author Martin R. Smith
#' @aliases StringToPhydat
#' @importFrom phangorn phyDat
#' @export
StringToPhyDat <- StringToPhydat <- function (string, tips, byTaxon = TRUE) {
  string <- strsplit(gsub('\\s+', '', string), '')[[1]]
  string <- matrix(string[string != '\n'], nrow=length(tips), byrow=byTaxon)
  rownames(string) <- tips
  phy <- phyDat(string, levels=c(which(as.character(0:9) %in% string) - 1, '-'), type='USER')
  phy
}

#' Extract character data from a phyDat object as a string
#' 
#' 
#' @param phy An object of class \code{\link{phyDat}}
#' @param ps Character specifying text, perhaps ';', to append to the end of the string
#' @param useIndex (default: TRUE) Print duplicate characters multiple 
#'        times, as they appeared in the original matrix
#' @param byTaxon If TRUE, write one taxon followed by the next.
#'                If FALSE, write one character followed by the next.
#' @param concatenate Logical specifying whether to concatenate all characters/taxa
#'                    into a single string, or to return a separate string
#'                    for each entry.
#' 
#' @author Martin R. Smith
#' @importFrom phangorn phyDat
#' @export
PhyToString <- function (phy, ps='', useIndex=TRUE, byTaxon=TRUE, concatenate=TRUE) {
  at <- attributes(phy)
  phyLevels <- at$allLevels		
  phyChars <- at$nr		
  phyContrast <- at$contrast == 1		
  phyIndex <- if (useIndex) at$index else seq_len(phyChars)		
  outLevels <- seq_len(ncol(phyContrast)) - 1		
  if (any(inappLevel <- phyLevels == '-')) outLevels[which(phyContrast[inappLevel])] <- '-'		

  levelLengths <- vapply(outLevels, nchar, integer(1))
  longLevels <- levelLengths > 1
  if (any(longLevels)) {
    if ('10' %in% outLevels && !(0 %in% outLevels)) {
      outLevels[outLevels == '10'] <- '0'
      longLevels['10'] <- FALSE
    }
    outLevels[longLevels] <- LETTERS[seq_len(sum(longLevels))]
  }
  
  levelTranslation <- apply(phyContrast, 1, function (x)		
    ifelse(sum(x) == 1, as.character(outLevels[x]), paste0(c('{', outLevels[x], '}'), collapse=''))		
  )		
  if (any(ambigToken <- apply(phyContrast, 1, all))) levelTranslation[ambigToken] <- '?'		
  ret <- vapply(phy, function (x) levelTranslation[x[phyIndex]], character(length(phyIndex)))		
  if (!byTaxon) ret <- t(ret) # Make each row correspond to a taxon
  ret <- if (concatenate || is.null(dim(ret))) { # If only one row, don't need to apply
    paste0(c(ret, ps), collapse='')		
  } else {		
    paste0(apply(ret, 1, paste0, collapse=''), ps)		
  }
  # Return:
  ret
}

#' @name AsBinary
#' @aliases AsBinary
#' @title Convert a number to binary
#' @description Provides a (reversed) binary representation of a decimal integer
#' @usage AsBinary(x)
#'
#' @param x Decimal integer to be converted to binary bits
#' 
#' @details 
#' Provides an array corresponding to binary digits 1, 2, 4, 8, 16, ...
#' 
#' Binary number 0100 (= decimal 4) will be represented as 0 0 1.
#' 
#' @return 
#' An array corresponding to binary digits 1, 2, 4, 8, 16, ...
#' 
#' 'Leading zeros' are not included.
#' 
#' @author Martin R. Smith, adapted from code posted to R mailing list by Spencer Graves
#' 
#' @examples
#'   AsBinary(4)  # 0 0 1
#'   AsBinary(10) # 0 1 0 1
#' 
#' @export
AsBinary <- function(x) {
	N <- length(x)
	xMax <- max(x)	
	ndigits <- (floor(logb(xMax, base=2))+1)
	Base.b <- array(NA, dim=c(N, ndigits))
	for (i in 1:ndigits){#i <- 1
		Base.b[, i] <- (x %% 2)
		x <- (x %/% 2)
	}
	if(N == 1) Base.b[1, ] else Base.b
}