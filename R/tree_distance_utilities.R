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

#' List matched clades as text
#' @param clades1,clades2 Logical matrices with columns specifying membership
#' of each corresponding matched clade 
#' @param taxonNames Character vector listing names of taxa corresponding to
#'  each row in `clades#`
#' @seealso VisualizeMatching
#' @keywords internal
#' @export
ReportMatching <- function (clades1, clades2, taxonNames) {
  clades1 <- apply(clades1, 2, function (x) paste0(
    paste(taxonNames[x], collapse=' '), ' : ', 
    paste(taxonNames[!x], collapse=' ')))
  clades2 <- apply(clades2, 2, function (x) paste0(
    paste(taxonNames[x], collapse=' '), ' : ', 
    paste(taxonNames[!x], collapse=' ')))
  
  # Return:
  paste(clades1, '=>', clades2)
}

#' Visualise a matching
#' 
#' Depicts the bipartitions that are matched between two trees using a 
#' specified Generalized Robinson Foulds tree distance measure.
#' 
#' @param Func Function used to construct tree similarity.
#' @param tree1,tree2 Trees of class `phylo`, with tips labelled identically.
#' @param setPar Logical specifying whether graphical parameters should be 
#' set to display trees side by side.
#' @param precision Integer specifying number of significant figures to display
#' when reporting matching scores.
#' @param Plot Function to use to plot trees.
#' @param \dots Additional parameters to send to `Plot`.
#' 
#' @author Martin R. Smith
#' @importFrom ape nodelabels edgelabels plot.phylo
#' @importFrom colorspace qualitative_hcl
#' @importFrom graphics par
#' @export
VisualizeMatching <- function(Func, tree1, tree2, setPar = TRUE, 
                              precision=3, Plot = plot.phylo, ...) {
  
  splits1 <- Tree2Splits(tree1)
  edge1 <- tree1$edge
  child1 <- edge1[, 2]
  partitionEdges1 <- vapply(colnames(splits1), 
                            function (node) which(child1 == node), integer(1))
  
  splits2 <- Tree2Splits(tree2)
  edge2 <- tree2$edge
  child2 <- edge2[, 2]
  partitionEdges2 <- vapply(colnames(splits2), 
                            function (node) which(child2 == node), integer(1))
  
  matching <- Func(tree1, tree2, reportMatching = TRUE)
  pairings <- attr(matching, 'matching')
  scores <- attr(matching, 'pairScores')
  pairScores <- signif(mapply(function (i, j) scores[i, j],
                              seq_along(pairings), pairings), precision)
  
  nSplits <- length(pairings)
  palette <- qualitative_hcl(nSplits, c=42, l=88)
  adjNo <- c(0.5, -0.2)
  adjVal <- c(0.5, 1.1)
  
  if (setPar) origPar <- par(mfrow=c(2, 1), mar=rep(0.5, 4))
  
  Plot(tree1)#, ...)
  edgelabels(seq_along(pairings), partitionEdges1[pairings], bg=palette, adj=adjNo)
  #edgelabels(seq_along(pairings), partitionEdges1, cex=0.8, font=2, frame='n', adj=adjVal)
  edgelabels(pairScores, partitionEdges1[pairings], frame='n', adj=adjVal, cex=0.8)
  
  Plot(tree2)#, ...)
  edgelabels(seq_along(pairings), partitionEdges2, bg=palette, adj=adjNo)
  #edgelabels(seq_along(pairings), partitionEdges2, cex=0.8, font=2, frame='n', adj=adjVal)
  edgelabels(pairScores, partitionEdges2, frame='n', adj=adjVal, cex=0.8)
  
  if (setPar) par(origPar)
}
