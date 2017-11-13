#' Search for most parsimonious trees
#'
#' Run standard search algorithms (\acronym{NNI}, \acronym{SPR} or \acronym{TBR}) 
#' to search for a more parsimonious tree.
#'  
#' @param tree a fully-resolved starting tree in \code{\link{phylo}} format, with the desired outgroup; edge lengths are not supported and will be deleted;
#' @template datasetParam
#' @template EdgeSwapperParam
#' @param maxIter the maximum number of iterations to perform before abandoning the search.
#' @param maxHits the maximum times to hit the best pscore before abandoning the search.
#' @param forestSize the maximum number of trees to return - useful in concert with \code{\link{consensus}}.
#'
#' @template InitializeDataParam
#' @template CleanUpDataParam
#' @template TreeScorerParam
#'
#' @template verbosityParam
#' @template treeScorerDots
#' 
#' @return {
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
#' \item \code{\link{Sectorial}}, alternative heuristic, useful for larger trees;
#' \item \code{\link{Ratchet}}, alternative heuristic, useful to escape local optima.
#' }
#'
#' @examples
#' data('Lobo')
#' njtree <- NJTree(Lobo.phy)
#'
#' \dontrun{
#' TreeSearch(njtree, Lobo.phy, maxIter=20, Rearrange=NNI)
#' TreeSearch(njtree, Lobo.phy, maxIter=20, Rearrange=SPR)
#' TreeSearch(njtree, Lobo.phy, maxIter=20, Rearrange=TBR)}
#' 
#' @keywords  tree 
#' 
#' @export
TreeSearch <- function (tree, dataset,
                        InitializeData = PhyDat2Morphy,
                        CleanUpData    = UnloadMorphy,
                        TreeScorer     = MorphyLength,
                        EdgeSwapper    = RootedTBRSwap,
                        maxIter = 100, maxHits = 20, forestSize = 1,
                        verbosity = 1, ...) {
  epsilon <- 1e-07
  hits <- 0L
  # initialize tree and data
  if (dim(tree$edge)[1] != 2 * tree$Nnode) stop("tree must be bifurcating; try rooting with ape::root")
  tree <- RenumberTips(tree, names(dataset))
  edgeList <- MatrixToList(tree$edge)
  edgeList <- RenumberEdges(edgeList[[1]], edgeList[[2]])

  initializedData <- InitializeData(dataset)
  on.exit(initializedData <- CleanUpData(initializedData))

  if (exists("forestSize") && length(forestSize) && forestSize > 1) {
    # TODO ForestSize > 1 not supported
    forest <- empty.forest <- vector('list', forestSize)
    forest[[1]] <- edgeList
  } else {
    forestSize <- 1 
  }
  bestScore <- if (is.null(attr(tree, 'score'))) {
    TreeScorer(edgeList[[1]], edgeList[[2]], initializedData, ...)
  } else {
    attr(tree, 'score')
  }
  if (verbosity > 0) cat("\n  - Performing tree search.  Initial score:", bestScore)
  returnSingle <- !(forestSize > 1)
  
  for (iter in 1:maxIter) {
    candidateEdges <- RearrangeEdges(edgeList[[1]], edgeList[[2]], dataset=initializedData, 
                             TreeScorer=TreeScorer, hits=hits, inputScore=bestScore,
                             EdgeSwapper=EdgeSwapper, minScore=bestScore,
                             returnSingle=returnSingle, iter=iter,
                             verbosity=verbosity, ...)
    iterScore <- candidateEdges[[3]]
    if (length(forestSize) && forestSize > 1) {
      ### TODO
      ###hits <- attr(trees, 'hits')
      ###if (iterScore == bestScore) {
      ###  forest[(hits - length(trees) + 1):hits] <- trees
      ###  tree <- sample(forest[1:hits], 1)[[1]]
      ###  attr(tree, 'score') <- iterScore
      ###  attr(tree, 'hits') <- hits
      ###} else if (iterScore < bestScore) {
      ###  bestScore <- iterScore
      ###  forest <- empty.forest
      ###  forest[seq_len(hits)] <- trees
      ###  tree <- sample(trees, 1)[[1]]
      ###  attr(tree, 'score') <- iterScore
      ###  attr(tree, 'hits') <- hits
      ###}      
    } else {
      if (iterScore < bestScore + epsilon) {
        hits <- candidateEdges[[4]]
        bestScore <- iterScore
        edgeList <- candidateEdges
      }
    }
    if (hits >= maxHits) break
  }
  if (verbosity > 0) cat("\n  - Final score", bestScore, "found", hits, "times after", iter, "iterations\n")  
  if (forestSize > 1) {
    if (hits < forestSize) forest <- forest[-((hits+1):forestSize)]
    attr(forest, 'hits') <- hits
    attr(forest, 'score') <- bestScore
    return(unique(forest))
  } else {
    tree$edge <- ListToMatrix(edgeList)
    attr(tree, 'hits') <- hits
    attr(tree, 'score') <- bestScore
    return(tree)
  }
}