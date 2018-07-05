#' Frequency of splits
#' @param referenceTree A tree of class phylo
#' @param forest a list of trees of class phylo, or a multiPhylo object
#' @return Number of trees in `forest` that contain each split in `referenceTree`
#' 
#' @author Martin R. Smith
#' @importFrom phangorn Descendants
#' @export
SplitFrequency <- function(referenceTree, forest) {
  tipIndex <-referenceTree$tip.label
  nTip <- length(tipIndex)
  powersOf2 <- 2 ^ (seq_len(nTip) - 1)
  if (class(forest) == 'phylo') forest <- list(forest)
  
  SplitNumber <- function (tips, tr) {
    included <- tipIndex %in% tr$tip.label[tips]
    as.character(min(c(sum(powersOf2[included]), sum(powersOf2[!included]))))
  }
  
  treeSplits <- Descendants(referenceTree, nTip + seq_len(Nnode(referenceTree)),
                            type='tips')
  splitNumbers <- vapply(treeSplits, SplitNumber, character(1), referenceTree)
  
  forestSplits <- table(vapply(forest, function (tr) {
    vapply(Descendants(tr, nTip + seq_len(nTip - 1L), type='tips'),
           SplitNumber, character(1), tr)
  }, character(nTip - 1L)))
  
  occurrences <- forestSplits[splitNumbers]
  # Root split appears twice!
  occurrences[duplicated(splitNumbers)] <- occurrences[2] <- occurrences[2] / 2L
  occurrences
}

#' Support colour
#' @param support A vector of doubles in the range 0-1
#' @param show1 Logical specifying whether to display values of 1. 
#'              A transparent white will be returned if `FALSE`.  
#' @return A string containing the hexadecimal code for a colour picked from a
#'         diverging scale, or `red` if a value is invalid.
#' @importFrom colorspace diverge_hcl
#' @export
SupportColour <- function (support, show1=TRUE) {
  # continuousScale <- rev(colorspace::heat_hcl(101, h=c(300, 75), c.=c(35, 95), l=c(15, 90), power=c(0.8, 1.2))) # Viridis prefered
  divergingScale <- rev(diverge_hcl(101, h=c(260, 0), c=100, l=c(50, 90), power=1.0))
  ifelse(is.na(support) | support < 0 | support > 1 | support == '', 'red',
         ifelse(support == 1 & !show1, "#ffffff00", divergingScale[(support * 100) + 1L]))
}

#' @export
SupportColor <- SupportColour

