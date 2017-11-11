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
#' @param dataset data table, perhaps from read.nexus.data
#' @param levels tokens - values that all characters migkt take
#' @param compress Compress identical transformation series into a single row of the phyDat object
#' For simplicity I have not retained support for contrast matrices or ambiguity.
#' @return a \code{phyDat} object
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
#' @return a dataset of class 'profileDat'
#'
#' @author Martin R. Smith; written with reference to phangorn:::prepareDataFitch
#' @export
PrepareDataProfile <- function (dataset, precision = 1e+06, warn = TRUE) {
  at <- attributes(dataset)
  nam <- at$names
  nLevel <- length(at$level)
  nChar <- at$nr
  cont <- attr(dataset, "contrast")
  nTip <- length(dataset)
  
  at$names <- NULL
  powers.of.2 <- 2L ^ c(0L:(nLevel - 1L))
  tmp <- cont %*% powers.of.2
  tmp <- as.integer(tmp)
  dataset <- unlist(dataset, recursive=FALSE, use.names=FALSE)
  ret <- tmp[dataset]
  ret <- as.integer(ret)
  attributes(ret) <- at
  inappLevel <- which(at$levels == "-")
  attr(ret, 'inappLevel') <- 2 ^ (inappLevel - 1)
  attr(ret, 'dim') <- c(nChar, nTip)  
  applicableTokens <- setdiff(powers.of.2, 2 ^ (inappLevel - 1))
  attr(ret, 'split.sizes') <- apply(ret, 1, function(x) vapply(applicableTokens, function (y) sum(x == y), integer(1)))
  attr(ret, 'info.amounts') <- InfoAmounts(ret, precision, warn=warn)
  attr(ret, 'bootstrap') <- c('info.amounts', 'split.sizes')
  dimnames(ret) <- list(NULL, nam)
  class(ret) <- 'profileDat'
  ret
}
