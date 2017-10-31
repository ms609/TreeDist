#' Search for most parsimonious trees
#'
#' Run standard search algorithms (\acronym{NNI}, \acronym{SPR} or \acronym{TBR}) 
#' to search for a more parsimonious tree.
#'  
#' @param tree a fully-resolved starting tree in \code{\link{phylo}} format, with the desired outgroup; edge lengths are not supported and will be deleted;
#' @template datasetTreeScorerParams
#' @param outgroup a vector listing the taxa in the outgroup;
#' @param concavity concavity constant for implied weighting (not currently implemented!); 
#' @param Rearrange rearrangement function to use; perhaps one of \kbd{RootedNNI}, \kbd{SPR}, or \kbd{TBR};
#' @param maxIter the maximum number of iterations to perform before abandoning the search;
#' @param maxHits the maximum times to hit the best pscore before abandoning the search;
#' @param forestSize the maximum number of trees to return - useful in concert with \code{\link{consensus}};
#' @template verbosityParam
#' @template treeScorerDots
#' 
#' @return{
#' This function returns a tree, with an attribute \code{pscore} conveying its parsimony score.
#' Note that the parsimony score will be inherited from the tree's attributes, which is only valid if it 
#' was generated using the same \code{data} that is passed here.
#' }
#' @author Martin R. Smith
#'
#' @seealso
#' \itemize{
#' \item \code{\link{Fitch}}, calculates parsimony score;
#' \item \code{\link{RootedNNI}}, conducts tree rearrangements;
#' \item \code{\link{SectorialSearch}}, alternative heuristic, useful for larger trees;
#' \item \code{\link{Ratchet}}, alternative heuristic, useful to escape local optima.
#' }
#'
#' @examples
#' data('Lobo')
#' outgroup <- c('Lingula', 'Mickwitzia', 'Neocrania')
#' njtree <- ape::root(nj(dist.hamming(Lobo.phy)), outgroup, resolve.root=TRUE)
#' njtree$edge.length <- NULL; njtree<-SetOutgroup(njtree, outgroup)
#'
#' \dontrun{
#' TreeSearch(njtree, Lobo.phy, outgroup, maxIter=20, Rearrange=NNI)
#' TreeSearch(njtree, Lobo.phy, outgroup, maxIter=20, Rearrange=SPR)
#' TreeSearch(njtree, Lobo.phy, outgroup, maxIter=20, Rearrange=TBR)}
#' 
#' @keywords  tree 
#' 
#' @export
TreeSearch <- function (tree, dataset, 
                        InitializeData = InitFitch,
                        CleanUpData = DestroyFitch,
                        TreeScorer = IFitchScore,
                        Rearrange = RootedTBR,
                        maxIter = 100, maxHits = 20, forestSize = 1,
                        verbosity = 1, ...) {
  # initialize tree and data
  InitializeData(tree, dataset)
  on.exit(CleanUpData)
  if (is.null(treeOrder <- attr(tree, 'order')) || treeOrder != 'preorder') tree <- Preorder(tree)
  tree$edge.length <- NULL # Edge lengths are not supported
  attr(tree, 'hits') <- 1
  if (exists("forestSize") && length(forestSize) && forestSize > 1) {
    forest <- empty.forest <- vector('list', forestSize)
    forest[[1]] <- tree
  } else {
    forestSize <- 1 
  }
  if (is.null(attr(tree, 'score'))) attr(tree, 'score') <- TreeScorer(tree=tree, ...)
  bestScore <- attr(tree, 'score')
  if (verbosity > 0) cat("\n  - Performing tree search.  Initial score:", bestScore)
  returnSingle <- !(forestSize > 1)
  
  for (iter in 1:maxIter) {
    trees <- RearrangeTree(tree, TreeScorer, Rearrange, minScore=bestScore,
                           returnSingle=returnSingle, iter=iter, verbosity=verbosity, ...)
    iterScore <- attr(trees, 'score')
    if (length(forestSize) && forestSize > 1) {
      hits <- attr(trees, 'hits')
      if (iterScore == bestScore) {
        forest[(hits - length(trees) + 1):hits] <- trees
        tree <- sample(forest[1:hits], 1)[[1]]
        attr(tree, 'score') <- iterScore
        attr(tree, 'hits') <- hits
      } else if (iterScore < bestScore) {
        bestScore <- iterScore
        forest <- empty.forest
        forest[seq_len(hits)] <- trees
        tree <- sample(trees, 1)[[1]]
        attr(tree, 'score') <- iterScore
        attr(tree, 'hits') <- hits
      }      
    } else {
      if (iterScore <= bestScore) {
        bestScore <- iterScore
        tree <- trees
      }
    }
    if (attr(trees, 'hits') >= maxHits) break
  }
  if (verbosity > 0) cat("\n  - Final score", attr(tree, 'score'), "found", attr(tree, 'hits'), "times after", iter, "iterations\n")  
  if (forestSize > 1) {
    if (hits < forestSize) forest <- forest[-((hits+1):forestSize)]
    attr(forest, 'hits') <- hits
    attr(forest, 'score') <- bestScore
    return(unique(forest))
  } else {
    return(tree)
  }
}