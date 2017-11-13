#' @describeIn TreeSearch Tree Search from Edge lists
#' @template edgeListParam
#' @template dataForFunction
#' @author Martin R. Smith
#' @keywords internal
#' @export
EdgeListSearch <- function (edgeList, dataset,
                          TreeScorer = MorphyLength,
                          EdgeSwapper = RootedTBRSwap,
                          maxIter=100, maxHits=20, 
                          bestScore=NULL,
                          stopAtScore=NULL, forestSize=1L, 
                          cluster=NULL, verbosity=1L, ...) {
  epsilon <- 1e-07                        
  if (!is.null(forestSize) && length(forestSize)) {
    if (forestSize > 1L) {
      stop("TODO: Forests not supported")
      #### forest <- empty.forest <- vector('list', forestSize)
      #### forest[[1]] <- edgeList
    } else {
      forestSize <- 1L
    }
  }
  if (is.null(bestScore)) {
    if (length(edgeList) < 3L) {
      bestScore <- TreeScorer(edgeList[[1]], edgeList[[2]], dataset)
    } else {
      bestScore <- edgeList[[3]]
    }
  }
  if (verbosity > 0L) cat("\n  - Performing tree search.  Initial score:", bestScore)
  if (!is.null(stopAtScore) && bestScore < stopAtScore + epsilon) return(edgeList)
  returnSingle <- !(forestSize > 1L)
  hits <- 0L
  
  for (iter in 1:maxIter) {
    candidateLists <- RearrangeEdges(edgeList[[1]], edgeList[[2]], dataset=dataset, 
                             TreeScorer=TreeScorer, hits=hits, inputScore=bestScore,
                             EdgeSwapper=EdgeSwapper, minScore=bestScore,
                             returnSingle=returnSingle, iter=iter,
                             verbosity=verbosity, ...)
    scoreThisIteration <- candidateLists[[3]]
    hits <- candidateLists[[4]]
    if (forestSize > 1L) {
      stop("TODO re-code this")
      ###if (scoreThisIteration == bestScore) {
      ###  forest[(hits-length(candidateLists)+1L):hits] <- candidateLists ## TODO Check that length still holds
      ###  edgeList  <- sample(forest[1:hits], 1)[[1]]
      ###  bestScore <- scoreThisIteration
      ###  hits      <- hits + 1L
      ###} else if (scoreThisIteration < bestScore) {
      ###  bestScore <- scoreThisIteration
      ###  forest <- empty.forest
      ###  forest[1:hits] <- candidateLists
      ###  edgeList <- sample(candidateLists , 1)[[1]]
      ###  attr(edgeList, 'score') <- scoreThisIteration
      ###}
    } else {
      if (scoreThisIteration < bestScore + epsilon) {
        bestScore <- scoreThisIteration
        edgeList  <- candidateLists
        if (!is.null(stopAtScore) && bestScore < stopAtScore + epsilon) return(edgeList)
      }
    }
    if (hits >= maxHits) break
  }
  if (verbosity > 0L) cat("\n  - Final score", bestScore, "found", hits, "times after", iter, "rearrangements\n")  
  if (forestSize > 1L) {
    if (hits < forestSize) forest <- forest[-((hits+1):forestSize)]
    attr(forest, 'hits') <- hits
    attr(forest, 'score') <- bestScore
    return (unique(forest))
  } else {
    edgeList[3:4] <- c(bestScore, hits)
    return(edgeList)
  }
}

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
#' TreeSearch(njtree, Lobo.phy, maxIter=20, EdgeSwapper=NNISwap)
#' TreeSearch(njtree, Lobo.phy, maxIter=20, EdgeSwapper=RootedSPRSwap)
#' TreeSearch(njtree, Lobo.phy, maxIter=20, EdgeSwapper=TBRSwap)}
#' 
#' @keywords  tree 
#' 
#' @export
TreeSearch <- function (tree, dataset,
                        InitializeData = PhyDat2Morphy,
                        CleanUpData    = UnloadMorphy,
                        TreeScorer     = MorphyLength,
                        EdgeSwapper    = RootedTBRSwap,
                        maxIter = 100L, maxHits = 20L, forestSize = 1L,
                        verbosity = 1L, ...) {
  # initialize tree and data
  if (dim(tree$edge)[1] != 2 * tree$Nnode) stop("tree must be bifurcating; try rooting with ape::root")
  tree <- RenumberTips(tree, names(dataset))
  edgeList <- MatrixToList(tree$edge)
  edgeList <- RenumberEdges(edgeList[[1]], edgeList[[2]])

  initializedData <- InitializeData(dataset)
  on.exit(initializedData <- CleanUpData(initializedData))

  bestScore <- attr(tree, 'score')
  edgeList <- EdgeListSearch(edgeList, initializedData, TreeScorer=TreeScorer, EdgeSwapper=EdgeSwapper,
                 maxIter = maxIter, maxHits = maxHits, forestSize = forestSize, verbosity = verbosity)
  
  tree$edge <- ListToMatrix(edgeList)
  attr(tree, 'score') <- edgeList[[3]]
  attr(tree, 'hits') <- edgeList[[4]]
  # Return:
  tree 
}