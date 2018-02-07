## Copied from phangorn file phyDat.R.  Edited for style only.
FastTable <- function (dataset) {                                                                                 
  if(!is.data.frame(dataset)) {
    dataset <- as.data.frame(dataset, stringsAsFactors = FALSE)                    
  }
  da <- do.call("paste", c(dataset, sep = "\r"))
  ind <- !duplicated(da)
  levels <- da[ind]
  cat <- factor(da, levels = levels)
  nl <- length(levels(cat))
  bin <- (as.integer(cat) - 1)                            
  pd <- nl
  bin <- bin[!is.na(bin)]
  if (length(bin)) bin <- bin + 1
  y <- tabulate(bin, pd)
  result <- list(index = bin, weights = y, dataset = dataset[ind,])
  result                                                                              
}   

#' Load phyDat object
#'
#' A convenient wrapper for \pkg{phangorn}'s \code{phyDat}
#'
#' @param dataset data table, perhaps from read.nexus.data.
#' @param levels tokens - values that all characters might take.
#' @param compress Compress identical transformation series into a single row of the phyDat object.
#' For simplicity I have not retained support for contrast matrices or ambiguity.
#'
#' @return a \code{phyDat} object.
#'
#' @importFrom stats na.omit
#' @export
PhyDat <- function (dataset, levels = NULL, compress = TRUE) {
  if (is.null(levels)) stop("Levels not supplied")
  nam <- names(dataset)
  # dataset <- as.data.frame(t(dataset), stringsAsFactors = FALSE)
  if (length(dataset[[1]]) == 1) {
      compress <- FALSE
  }
  if (compress) {
    ddd     <- FastTable(dataset)
    dataset <- ddd$dataset
    weight  <- ddd$weight
    index   <- ddd$index
    n.rows  <- length(dataset[[1]])
  } else {
    n.rows <- length(dataset[[1]])
    weight <- rep(1, n.rows)
    index <- 1:n.rows
  }
  n.levels <- length(levels)
  contrast <- diag(n.levels)
  all.levels <- levels
  att <- attributes(dataset)
  dataset <- lapply(dataset, match, all.levels)
  attributes(dataset) <- att
  row.names(dataset) <- as.character(1:n.rows)
  dataset <- na.omit(dataset)
  tmp  <- match(index, attr(dataset, "na.action"))
  index <- index[is.na(tmp)]
  index <- match(index, unique(index))
  rn <- as.numeric(rownames(dataset))
  attr(dataset, "na.action") <- NULL
  weight <- weight[rn]
  n.rows <- dim(dataset)[1]
  names(dataset) <- nam
  attr(dataset, "row.names") <- NULL
  attr(dataset, "weight")    <- weight
  attr(dataset, "nr")        <- n.rows
  attr(dataset, "nc")        <- length(levels)
  attr(dataset, "index")     <- index
  attr(dataset, "levels")    <- levels
  attr(dataset, "allLevels") <- all.levels
  attr(dataset, "type")      <- "USER"
  attr(dataset, "contrast")  <- contrast
  class(dataset)             <- "phyDat"
  dataset
}

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
  
  powers.of.2 <- 2L ^ c(0L:(nLevel - 1L))
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
  
  powers.of.2 <- 2L ^ c(0L:(nLevel - 1L))
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
#' Smallest number of steps that a character can take on any tree
#' 
#' @param states Integer vector listing the tokens that may be present at each tip along a single character,
#'               with each token represented as a binary digit; e.g. a value of 11 means that
#'               the tip may have tokens 0, 1 or 3 (as 11 = 2^0 + 2^1 + 2^3).
#'               As the minimum steps can be found when inapplicables occur together,
#'               inapplicable tokens can be denoted as ?s or with the integer 0 (not 2^0).
#'               
#' @return An integer specifying the minimum number of steps that the character must contain
MinimumSteps <- function (states) {
  
  tokens <- AsBinary(unique(states[states>0])) > 0
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
    
    tokens <- tokens[statesRemaining, !tokenNecessary]
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
