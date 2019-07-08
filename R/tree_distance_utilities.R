#' Wrapper for tree distance calculations
#' 
#' Calls tree distance functions from trees or lists of trees
#' 
#' @inheritParams MutualArborealInfo
#' @param Func Tree distance function.
#' @param \dots Additional arguments to `Func`.
#' 
#' @author Martin R. Smith
#' @keywords internal
#' @export
CalculateTreeDistance <- function (Func, tree1, tree2, reportMatching, ...) {
  if (class(tree1) == 'phylo') {
    if (class(tree2) == 'phylo') {
      if (length(setdiff(tree1$tip.label, tree2$tip.label)) > 0) {
        stop("Tree tips must bear identical labels")
      }
      
      Func(Tree2Splits(tree1), Tree2Splits(tree2), reportMatching, ...)
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

ReportMatching <- function (clades1, clades2, taxonNames) {
  clades1 <- apply(clades1, 2, function (x) paste0(
    paste(taxonNames[x], collapse=' '), ':', 
    paste(taxonNames[!x], collapse=' ')))
  clades2 <- apply(clades2, 2, function (x) paste0(
    paste(taxonNames[x], collapse=' '), ':', 
    paste(taxonNames[!x], collapse=' ')))
  
  # Return:
  paste(clades1, '=>', clades2)
}
