#' Wrapper for tree distance calculations
#' 
#' Calls tree distance functions from trees or lists of trees
#' 
#' @inheritParams MutualPhylogeneticInfo
#' @param Func Tree distance function.
#' @param \dots Additional arguments to `Func`.
#' 
#' @author Martin R. Smith
#' @keywords internal
#' @importFrom TreeTools as.Splits
#' @export
CalculateTreeDistance <- function (Func, tree1, tree2, reportMatching, ...) {
  if (class(tree1) == 'phylo') {
    labels1 <- tree1$tip.label
    if (class(tree2) == 'phylo') {
      if (length(setdiff(labels1, tree2$tip.label)) > 0) {
        stop("Tree tips must bear identical labels")
      }
      
      Func(as.Splits(tree1, asSplits = FALSE),
           as.Splits(tree2, tipLabels = labels1, asSplits = FALSE),
           reportMatching = reportMatching, ...)
    } else {
      splits1 <- as.Splits(tree1, tipLabels = labels1, asSplits = FALSE)
      vapply(tree2,
             function (tr2) Func(splits1, 
                                 as.Splits(tr2, tipLabels = labels1, 
                                           asSplits = FALSE), ...),
             double(1))
    }
  } else {
    if (class(tree2) == 'phylo') {
      splits1 <- as.Splits(tree2, asSplits = FALSE)
      labels1 <- tree2$tip.label
      vapply(tree1,
             function (tr2) Func(splits1,
                                 as.Splits(tr2, tipLabels = labels1,
                                           asSplits = FALSE), ...),
             double(1))
    } else {
      splits1 <- as.Splits(tree1, asSplits = FALSE)
      splits2 <- as.Splits(tree2, tipLabels = tree1[[1]]$tip.label,
                           asSplits = FALSE)
      matrix(mapply(Func, rep(splits2, each=length(splits1)), splits1, ...), 
             length(splits1), length(splits2),
             dimnames = list(names(tree1), names(tree2)))
    }
  }
}

#' Entropy in bits
#' 
#' Reports the entropy of a vector of probabilities, in bits.
#' Probabilities should sum to one.  Probabilities equalling zero will be 
#' ignored.
#' 
#' @param p Numeric vector specifying probabilities of outcomes.
#' 
#' @examples
#' Entropy(rep(0.5, 2)) # = 1
#' Entropy(c(1/4, 1/4, 0, 1/4, 1/4)) # = 2
#' 
#' @return Entropy of the specified probabilities, in bits
#' @author Martin R. Smith
#' @export
Entropy <- function (p) -sum(p[p > 0] * log2(p[p > 0]))


#' Normalize information against total present in both starting trees
#' @param unnormalized Numeric value to be normalized.
#' @param tree1,tree2 Trees from which `unnormalized` was calculated
#' @param InfoInTree Function to calculate the information content of each tree
#' @param infoInBoth Numeric speecifying information content of both trees
#' independently (optional)
#' @param how Method for normalization
#' @param Func Function that takes as inputs `tree1Info` and `tree2Info`, and
#' returns a normalizing constant against which to divide `unnormalized`.
#' @param \dots Additional parameters to `InfoInTree`` or `how`.
#' @keywords internal
#' @author Martin R. Smith
#' @export
NormalizeInfo <- function (unnormalized, tree1, tree2, InfoInTree, 
                           infoInBoth = NULL,
                           how = TRUE, Combine = '+', ...) {
  
  CombineInfo <- function (tree1Info, tree2Info, Combiner = Combine) {
    if (length(tree1Info) == 1 || length(tree2Info) == 1) {
      mapply(Combiner, tree1Info, tree2Info)
    } else {
      outer(tree1Info, tree2Info, Combiner)
    }
  }
  
  if (mode(how) == 'logical') {
    if (how == FALSE) {
      return (unnormalized)
    } else {
      if (is.null(infoInBoth)) 
        infoInBoth <- CombineInfo(InfoInTree(tree1, ...), InfoInTree(tree2, ...))
    }
  } else if (mode(how) == 'function') {
    if (is.null(infoInBoth)) 
      infoInBoth <- CombineInfo(InfoInTree(tree1, ...), InfoInTree(tree2, ...),
                                Combiner = how)
  } else {
    infoInBoth <- how
  }
  
  # Return:
  unnormalized / infoInBoth
}

#' List clades as text
#' @param splits,splits1,splits2 Logical matrices with columns specifying membership
#' of each corresponding matched clade.
#' @return `ReportMatching` returns a character vector describing each pairing 
#' in a matching.
#'   
#' @seealso VisualizeMatching
#' @author Martin R. Smith
#' @keywords internal
#' @export
ReportMatching <- function (splits1, splits2, realMatch = TRUE) {
  paste(as.character(splits1), ifelse(realMatch, '=>', '..'), 
        as.character(splits2))
}
