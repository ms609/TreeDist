#' Prepare data for Profile Parsimony
#' 
#' @param dataset dataset of class \code{phyDat}
#' @param precision number of random trees to generate when calculating Profile curves. 
#'                  With 22 tokens (taxa):
#'                  - Increasing precision from 4e+05 to 4e+06 reduces error by a mean of 
#'                  0.005 bits for each step after the first (max = 0.11 bits, sd=0.017 bits)
#'                  - Increasing precision from 1e+06 to 4e+06 reduces error by a mean of 
#'                  0.0003 bits for each step after the first (max = 0.046 bits, sd=0.01 bits)
#'                  
#' @template warnParam
#'
#' @return An object of class phyDat with additional attributes:
#'         \code{info.amounts}: details the information represented by each character when subject to N 
#'         additional steps.
#'         \code{split.sizes}: The size of the splits implied by each character
#'         \code{bootstrap}: The character vector \code{c('info.amounts', 'split.sizes')}, indicating 
#'                           attributes to sample when bootstrapping the dataset 9e.g. in Ratchet searches).
#'
#' @author Martin R. Smith; written with reference to phangorn:::prepareDataFitch
#' @export
PrepareDataProfile <- function (dataset, precision = 1e+06, warn = TRUE) {
  at <- attributes(dataset)
  nLevel <- length(at$level)
  nChar <- at$nr
  cont <- attr(dataset, "contrast")
  nTip <- length(dataset)
  
  powers.of.2 <- 2L ^ c(0L, seq_len(nLevel - 1L))
  tmp <- cont %*% powers.of.2
  tmp <- as.integer(tmp)
  unlisted <- unlist(dataset, recursive=FALSE, use.names=FALSE)
  binaryMatrix <- tmp[unlisted]
  attr(binaryMatrix, 'dim') <- c(nChar, nTip)
  
  attr(dataset, 'info.amounts') <- InfoAmounts(binaryMatrix, precision, warn=warn)
  if (!any(attr(dataset, 'bootstrap') == 'info.amounts')) {
    attr(dataset, 'bootstrap') <- c(attr(dataset, 'bootstrap'), 'info.amounts')
  }
  
  ####  inappLevel <- which(at$levels == "-")
  ####  applicableTokens <- setdiff(powers.of.2, 2 ^ (inappLevel - 1))
  ####  
  ####  attr(dataset, 'split.sizes') <- apply(binaryMatrix, 1, function(x) {
  ####      vapply(applicableTokens, function (y) sum(x == y), integer(1))
  ####    })
  
  dataset
}


#' @describeIn PrepareDataProfile Prepare data for implied weighting
#' @export
PrepareDataIW <- function (dataset) {
  at <- attributes(dataset)
  nLevel <- length(at$level)
  nChar <- at$nr
  cont <- attr(dataset, "contrast")
  nTip <- length(dataset)
  
  powers.of.2 <- 2L ^ c(0L, seq_len(nLevel - 1L))
  inappLevel <- which(at$levels == "-")
  cont[, inappLevel] <- 0
  tmp <- as.integer(cont %*% powers.of.2)
  unlisted <- unlist(dataset, use.names=FALSE)
  binaryMatrix <- matrix(tmp[unlisted], nChar, nTip, byrow=FALSE)
  
  attr(dataset, 'min.steps') <- apply(binaryMatrix, 1, MinimumSteps)
  
  # Return:
  dataset
}

#' Minimum steps
#' 
#' The smallest number of steps that a character can take on any tree.
#' 
#' 
#' @param states Integer vector listing the tokens that may be present at each 
#' tip along a single character, with each token represented as a binary digit;
#' e.g. a value of 11 means that
#' the tip may have tokens 0, 1 or 3 (as 11 = 2^0 + 2^1 + 2^3).
#' As the minimum steps can be found when inapplicables occur together,
#' inapplicable tokens can be denoted as ?s or with the integer 0 (not 2^0).
#' Tokens that are ambiguous for an inapplicable and an applicable
#' state are not presently supported.
#'               
#' @return An integer specifying the minimum number of steps that the character
#'  must contain.
#'
#' @examples {
#'   data('inapplicable.datasets')
#'   myPhyDat <- inapplicable.phyData[[4]] # or as.phyDat(read.nexus.data('filepath'))
#'   class(myPhyDat) # phyDat object
#'   
#'   # Convert list of character codings to an array
#'   myData <- vapply(myPhyDat, I, myPhyDat[[1]])
#'   
#'   # Convert phyDat's representation of states to binary
#'   myContrast <- attr(myPhyDat, 'contrast')
#'   tokens <- colnames(myContrast)
#'   binaryContrast <- integer(length(tokens))
#'   tokenApplicable <- tokens != '-'
#'   binaryContrast[tokenApplicable] <- 2 ^ (seq_len(sum(tokenApplicable)) - 1)
#'   binaryValues <- apply(myContrast, 1, 
#'     function (row) sum(binaryContrast[as.logical(row)]))
#'   myStates <- matrix(binaryValues[myData], nrow=nrow(myData),
#'                      ncol=ncol(myData), dimnames=dimnames(myData))
#'  
#'   # Finally, work out minimum steps 
#'   apply(myStates, 1, MinimumSteps)
#'   
#' }
#'
#' @author Martin R. Smith
#' @export
MinimumSteps <- function (states) {
  
  uniqueStates <- unique(states[states>0])
  if (length(uniqueStates) < 2) return (0)
  tokens <- AsBinary(uniqueStates) > 0
  lastDim <- dim(tokens)
  tokensUsed <- 0
  
  repeat {
    tokens <- tokens[, !duplicated(t(tokens)), drop=FALSE]
    unambiguous <- rowSums(tokens) == 1
    tokenNecessary <- apply(tokens[unambiguous, , drop=FALSE], 2, any)
    statesRemaining <- !unambiguous
    statesRemaining[statesRemaining] <- rowSums(tokens[statesRemaining, tokenNecessary, drop=FALSE]) == 0
    tokensUsed <- tokensUsed + sum(tokenNecessary)
    
    if (!any(statesRemaining)) return (tokensUsed - 1)
    
    tokens <- tokens[statesRemaining, !tokenNecessary, drop=FALSE]
    if (identical(dim(tokens), lastDim)) {
      unnecessary <- colSums(tokens) == 1
      if (any(unnecessary)) {
        tokens <- tokens[, !unnecessary, drop=FALSE]
      } else {
        stop("The token configuration [", paste(states, collapse=" "), 
             "] is not correctly handled by MinimumSteps.\n Please report this bug at ",
             "https://github.com/ms609/TreeSearch/issues/new")
      }
    }
    lastDim <- dim(tokens)
  }
  
}
