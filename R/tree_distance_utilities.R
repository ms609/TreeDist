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
#' @importFrom TreeTools as.Splits TipLabels
#' @export
CalculateTreeDistance <- function (Func, tree1, tree2,
                                   reportMatching = FALSE, ...) {
  class1 <- class(tree1)
  single1 <- (class1 == 'phylo' || class1 == 'Splits')
  labels1 <- TipLabels(tree1)
  nTip <- length(labels1)
    
  class2 <- class(tree2)
  single2 <- (class2 == 'phylo' || class2 == 'Splits')
  labels2 <- TipLabels(tree2)
  
  if (length(setdiff(labels1, labels2)) > 0) {
    stop("Tree tips must bear identical labels")
  }
  
  if (single1) {
    if (single2) {
      .TreeDistanceOneOne(Func, tree1, tree2, tipLabels = labels1, nTip,
                          reportMatching = reportMatching)
    } else {
      .TreeDistanceOneMany(Func, oneSplit = tree1, manySplits = tree2,
                           tipLabels = labels1, nTip = nTip)
    }
  } else {
    if (single2) {
      .TreeDistanceOneMany(Func, oneSplit = tree2, manySplits = tree1,
                           tipLabels = labels2)
    } else {
      .TreeDistanceManyMany(Func, tree1, tree2,
                            tipLabels = labels1)
    }
  }
}

.TreeDistanceOneOne <- function (Func, split1, split2, tipLabels, 
                                 nTip = length(tipLabels), reportMatching, ...) {
  Func(as.Splits(split1, asSplits = reportMatching),
       as.Splits(split2, tipLabels = tipLabels, asSplits = reportMatching),
       nTip = nTip, reportMatching = reportMatching, ...)
}

.TreeDistanceOneMany <- function (Func, oneSplit, manySplits, 
                                  tipLabels, nTip = length(tipLabels), ...) {
  s1 <- as.Splits(oneSplit, tipLabels = tipLabels, asSplits = FALSE)
  vapply(manySplits,
         function (s2) Func(s1, as.Splits(s2, tipLabels = tipLabels,
                                          asSplits = FALSE),
                            nTip = nTip, ...),
         double(1))
}

.TreeDistanceManyMany <- function (Func, splits1, splits2, 
                                   tipLabels, nTip = length(tipLabels), ...) {
  if (identical(splits1, splits2)) {
    splits <- as.Splits(splits1, tipLabels = tipLabels, asSplits = FALSE)
    nSplits <- length(splits)
    notLastSplit <- nSplits - 1L
    ret <- matrix(0, nSplits, nSplits)
    is <- matrix(seq_len(nSplits), nSplits, nSplits)
    ret[upper.tri(ret)] <- mapply(Func,
                                  splits[t(is)[lower.tri(is)]],
                                  splits[is[lower.tri(is)]],
                                  nTip = nTip,
                                  reportMatching = FALSE)
    ret[lower.tri(ret)] <- t(ret)[lower.tri(ret)]
    
    # Return:
    ret
  } else {
    splits1 <- as.Splits(splits1, tipLabels = tipLabels, asSplits = FALSE)
    splits2 <- as.Splits(splits2, tipLabels = tipLabels, asSplits = FALSE)
    matrix(mapply(Func, rep(splits2, each=length(splits1)), splits1,
                  nTip = nTip, ...),
           length(splits1), length(splits2),
           dimnames = list(names(splits1), names(splits2)))
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
