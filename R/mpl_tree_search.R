#' @describeIn Ratchet returns a list of optimal trees produced by nSearch Ratchet searches
#' @param nSearch Number of Ratchet searches to conduct (for RatchetConsensus)
#' @export
RatchetConsensus <- function (tree, dataset, ratchHits=10,
                              searchIter=500, searchHits=20, verbosity=0L, 
                              swappers=list(RootedNNISwap), nSearch=10, ...) {
  trees <- lapply(logical(nSearch), function (x) MorphyRatchet(tree, dataset, ratchIter=1, 
              searchIter=searchIter, searchHits=searchHits, verbosity=verbosity, swappers=swappers, ...))
  scores <- vapply(trees, function (x) attr(x, 'score'), double(1))
  trees <- unique(trees[scores == min(scores)])
  cat ("Found", length(trees), 'unique trees from ', nSearch, 'searches.')
  return (trees)
}


#' Tree Search
#'
#' Performs a tree search using a Morphy object
#' 
#' Does the hard work of searching for a most parsimonious tree, given the
#' parent and child vectors of a tree arranged in preorder (perhaps with 
#' \code{\link{RenumberEdges}}).
#' End-users are expected to access this function through its wrapper, TreeSearch
#' It is also called directly by MorphyRatchet and Sectorial functions
#'
#' @template edgeListParam
#' @template dataForFunction
#' @param EdgeSwapper Function to use to rearrange trees; example: 
#'                  \code{\link{TBRSwap}}.
#' @param maxIter maximum iterations to conduct.
#' @param maxHits stop search after this many hits.
#' @template stopAtScoreParam
#' @param forestSize how many trees to hold.
#' @template clusterParam
#' @template verbosityParam
#' @param \dots additional variables to pass to \code{\link{MorphyRearrange}}.
#'
#' @author Martin R. Smith
#' 
#' @keywords internal
#' @export

EdgeListSearch <- function (edgeList, dataset,
                          TreeScorer = MorphyLength,
                          EdgeSwapper = RootedTBRSwap,
                          #morphyObj, 
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
  hitsThisTime <- 0 
  if (verbosity > 0L) cat("\n  - Performing tree search.  Initial score:", bestScore)
  if (!is.null(stopAtScore) && bestScore < stopAtScore + epsilon) return(edgeList)
  returnSingle <- !(forestSize > 1L)
  
  for (iter in 1:maxIter) {
    candidateLists <- RearrangeEdges(edgeList[[1]], edgeList[[2]], dataset=dataset, 
                             TreeScorer=TreeScorer, hits=hits, inputScore=bestScore,
                             EdgeSwapper=EdgeSwapper, minScore=bestScore,
                             returnSingle=returnSingle, iter=iter,
                             verbosity=verbosity, ...)
    scoreThisIteration <- candidateLists[[3]]
    if (forestSize > 1L) {
      stop("TODO re-code this")
      ###if (scoreThisIteration == bestScore) {
      ###  forest[(hitsThisTime-length(candidateLists)+1L):hitsThisTime] <- candidateLists ## TODO Check that length still holds
      ###  edgeList  <- sample(forest[1:hitsThisTime], 1)[[1]]
      ###  bestScore <- scoreThisIteration
      ###  hitsThisTime      <- hitsThisTime + 1L
      ###} else if (scoreThisIteration < bestScore) {
      ###  bestScore <- scoreThisIteration
      ###  forest <- empty.forest
      ###  forest[1:hitsThisTime] <- candidateLists
      ###  edgeList <- sample(candidateLists , 1)[[1]]
      ###  attr(edgeList, 'score') <- scoreThisIteration
      ###}
    } else {
      if (scoreThisIteration < bestScore + epsilon) {
        hitsThisTime <- candidateLists[[4]] + if (bestScore < scoreThisIteration + epsilon) hitsThisTime else 0L
        bestScore <- scoreThisIteration
        edgeList  <- candidateLists
        if (!is.null(stopAtScore) && bestScore < stopAtScore + epsilon) return(edgeList)
      }
    }
    if (hitsThisTime >= maxHits) break
  }
  if (verbosity > 0L) cat("\n  - Final score", bestScore, "found", hitsThisTime, "times after", iter, "rearrangements\n")  
  if (forestSize > 1L) {
    if (hitsThisTime < forestSize) forest <- forest[-((hitsThisTime+1):forestSize)]
    attr(forest, 'hits') <- hitsThisTime
    attr(forest, 'score') <- bestScore
    return (unique(forest))
  } else {
    edgeList[3:4] <- c(bestScore, hitsThisTime)
    return(edgeList)
  }
}

###   #' Sectorial Search
###   #'
###   #' \code{SectorialSearch} performs a sectorial search on a tree, preserving the position of the root.
###   #'
###   #' \code{InapplicableSectorial} performs a sectorial search on the tree specified. A sectorial search 
###   #' detaches a random part of the tree, performs rearrangments on this subtree, then reattaches it 
###   #' to the main tree (Goloboff, 1999).
###   #' The improvement to local \var{score} hopefully (but not necessarily) improves the overall \var{score}.
###   #' As such, the output of \code{InapplicableSectorial} should be treated by further \acronym{TBR (/SPR/NNI)}
###   #' rearrangements and only retained if the ultimate parsimony score is better than 
###   #' that of the original tree.
###   #' 
###   #' \code{SectorialSearch} is a basic recipe that runs \code{InapplicableSectorial} followed by a few rounds
###   #' of tree rearrangement, returning a tree whose \var{score} is no worse than that of \code{start.tree}.
###   #' 
###   #' @param tree a rooted, resolved tree in \code{\link{phylo}} format from which to start the search;
###   #' @template datasetParam
###   #' @param maxIter maximum number of rearrangements to perform on each sectorial iteration;
###   #' @template verbosityParam
###   #' @param rearrangements method to use when rearranging subtrees: NNI, SPR or TBR;
###   #' @param \dots other arguments to pass to subsequent functions.
###   #' 
###   #' @return a rooted tree of class \code{phylo}.
###   #' 
###   #' @references Goloboff, P. (1999). \cite{Analyzing large data sets in reasonable times: solutions for composite optima.} Cladistics, 15(4), 415-428. doi:\href{http://dx.doi.org/10.1006/clad.1999.0122}{10.1006/clad.1999.0122}
###   #' 
###   #' @author Martin R. Smith
###   #' 
###   #' @seealso \code{\link{TreeSearch}}
###   #' @seealso \code{\link{MorphyRatchet}}
###   #' 
###   #' @examples
###   #' require('ape')
###   #' data('SigSut')
###   #' outgroup <- c('Lingula', 'Mickwitzia', 'Neocrania')
###   #' njtree <- ape::root(nj(dist.hamming(SigSut.phy)), outgroup, resolve.root=TRUE)
###   #' njtree$edge.length <- NULL; njtree<-ape::root(njtree, outgroup, resolve.root=TRUE)
###   #' InapplicableSectorial(njtree, SigSut.phy, ratchIter=1, maxIter=50, largest.sector=7)
###   #' \dontrun{SectorialSearch(njtree, SigSut.phy) # Will be time-consuming }
###   #' 
###   #' 
###   #' @keywords  tree 
###   #' @export
###   SectorialSearch <- function
###   (tree, dataset, SectorialRearrangements=NNI, maxIter=2000,
###    subsequentRearrangements = list(RootedNNI, RootedTBR, 
###     RootedSPR, RootedNNI), verbosity=3, ...) {
###     if (class(dataset) != 'phyDat') stop("dataset must be of class phyDat, not", class(dataset))
###     morphyObj <- PhyDat2Morphy(dataset)
###     on.exit(morphyObj <- UnloadMorphy(morphyObj))
###     tree <- RenumberTips(Renumber(tree), names(dataset))
###     if (is.null(attr(tree, "score"))) {
###       attr(tree, "score") <- MorphyTreeLength(tree, morphyObj, ...)
###     }
###     bestScore <- attr(tree, "score")
###     if (verbosity > 0) cat("* Initial score:", bestScore)
###   
###     if (class(subsequentRearrangements) == 'function') rearrangements <- list(rearrangements)
###     if (class(SectorialRearrangements) != 'function') stop("SectorialRearrangements must be a function, e.g. TreeSearch::NNI")
###     
###     bestScore <- attr(tree, 'score')
###     tree <- RenumberTips(Renumber(tree), names(dataset))
###     if (length(bestScore) == 0) bestScore <- Fitch(tree, dataset, ...)
###     sect <- MorphySectorial(tree, morphyObj, verbosity=verbosity-1, ratchIter=30, 
###       maxIter=maxIter, maxHits=15, smallest.sector=6, 
###       largest.sector=length(tree$edge[,2L])*0.25, rearrangements=rearrangements)
###     sect <- TreeSearch(sect, dataset, Rearrange=subsequentRearrangements, maxIter=maxIter, maxHits=30, cluster=cluster, verbosity=verbosity)
###     if (attr(sect, 'score') <= bestScore) {
###       return (sect)
###     } else return (tree)
###   }
###   