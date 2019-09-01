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
#' @importFrom TreeSearch Tree2Splits
#' @export
CalculateTreeDistance <- function (Func, tree1, tree2, reportMatching, ...) {
  if (class(tree1) == 'phylo') {
    if (class(tree2) == 'phylo') {
      if (length(setdiff(tree1$tip.label, tree2$tip.label)) > 0) {
        stop("Tree tips must bear identical labels")
      }
      
      Func(Tree2Splits(tree1), Tree2Splits(tree2), 
           reportMatching = reportMatching, ...)
    } else {
      splits1 <- Tree2Splits(tree1)
      vapply(tree2, 
             function (tr2) Func(splits1, Tree2Splits(tr2), ...),
             double(1))
    }
  } else {
    if (class(tree2) == 'phylo') {
      splits1 <- Tree2Splits(tree2)
      vapply(tree1, 
             function (tr2) Func(splits1, Tree2Splits(tr2), ...),
             double(1))
    } else {
      splits1 <- lapply(tree1, Tree2Splits)
      splits2 <- lapply(tree2, Tree2Splits)
      matrix(mapply(Func, rep(splits2, each=length(splits1)), splits1), 
             length(splits1), length(splits2),
             dimnames = list(names(tree1), names(tree2)), ...)
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
    if (how == FALSE) return (unnormalized)
    if (is.null(infoInBoth)) 
      infoInBoth <- CombineInfo(InfoInTree(tree1, ...), InfoInTree(tree2, ...))
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

#' Tree Distance Return
#' 
#' Generates the return value for Generalized Robinson-Foulds distances,
#' with (optionally) a report of the matching.
#' 
#' @param pairScores,reportMatching,swapSplits,taxonNames1,taxonNames2 Corresponding
#' parameters from within each function
#' @param Score function taking pairScores and matchedSplits to generate a 
#' matching score.
#' 
#' @keywords internal
#' @author Martin R. Smith
#' @export
TreeDistanceReturn <- function (pairScores, maximize = FALSE,
                                reportMatching, swapSplits, 
                                splits1, splits2,
                                taxonNames = NULL) {
  if (dim(pairScores)[1] == 1) {
    Optimal <- if (maximize) which.max else which.min
    optimalMatching <- structure(Optimal(pairScores), class='solve_LSAP')
  } else {
    optimalMatching <- solve_LSAP(pairScores, maximize)
  }
  
  ret <- sum(pairScores[matrix(c(seq_along(optimalMatching), optimalMatching),
                               ncol=2L)])
  
  if (reportMatching) {
    if (swapSplits) {
      optimalMatching <- structure(match(seq_len(dim(pairScores)[2]), optimalMatching),
                                   class='solve_LSAP')
      attr(ret, 'pairScores') <- t(pairScores)
    } else {
      attr(ret, 'pairScores') <- pairScores
    }
    
    if (!is.null(taxonNames)) {
      matchedSplits <- !is.na(optimalMatching)
      attr(ret, 'matchedSplits') <- 
        ReportMatching(splits1[, matchedSplits, drop=FALSE], 
                       splits2[, optimalMatching[matchedSplits], 
                               drop=FALSE], 
                       taxonNames)
    }
    attr(ret, 'matching') <- optimalMatching
  }
  # Return:
  ret
}

#' List clades as text
#' @param splits,splits1,splits2 Logical matrices with columns specifying membership
#' of each corresponding matched clade 
#' @param taxonNames Character vector listing names of taxa corresponding to
#'  each row in `splits`
#' @return {
#'   `ReportMatching` returns a character vector describing each pairing in a matching.
#'   
#'   `IdentifySplits` returns a character vector describing each split in a splits 
#'   object, named according to the column names in the splits object.
#' }
#' @seealso VisualizeMatching
#' @author Martin R. Smith
#' @keywords internal
#' @export
ReportMatching <- function (splits1, splits2, taxonNames) {
  paste(IdentifySplits(splits1, taxonNames), '=>', IdentifySplits(splits2, taxonNames))
}

#' @describeIn ReportMatching List the distribution of terminals represented by a single splits object.
#' @export
IdentifySplits <- function (splits, taxonNames = rownames(splits)) {
  apply(splits, 2, function (x) paste0(
    paste(taxonNames[x], collapse=' '), ' : ', 
    paste(taxonNames[!x], collapse=' ')))
}
