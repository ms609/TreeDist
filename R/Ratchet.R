#' Parsimony Ratchet
#'
#' \code{Ratchet} uses the parsimony ratchet (Nixon 1999) to search for a more parsimonious tree.
#'
#' @template treeParam 
#' @param dataset a dataset in the format required by TreeScorer.
#' @template InitializeDataParam
#' @template CleanUpDataParam
#' @template treeScorerParam
#' @param Bootstrapper Function to perform bootstrapped rearrangements of tree. 
#'                     First arguments will be an edgeList and a dataset, initialized using \code{InitializeData}
#'                     Should return a rearranged edgeList.
#' @template swappersParam
#' @param BootstrapSwapper Function such as \code{\link{RootedNNISwap}} to use to rearrange trees
#'                         within \code{Bootstrapper}.
#' @param returnAll Set to \code{TRUE} to report all MPTs encountered during the search, perhaps to analyze consensus.
#' @param ratchIter stop when this many ratchet iterations have been performed.
#' @param ratchHits stop when this many ratchet iterations have found the same best score.
#' @param searchIter maximum rearrangements to perform on each bootstrap or ratchet iteration.
#' @param searchHits maximum times to hit best score before terminating a tree search within a ratchet iteration.
#' @param bootstrapIter maximum rearrangements to perform on each bootstrap  iteration (default: \code{maxIter}).
#' @param bootstrapHits maximum times to hit best score on each bootstrap  iteration (default: \code{maxHits}).
#' @template stopAtScoreParam
#' @template verbosityParam
#' @param suboptimal retain trees that are suboptimal by this score. Defaults to 1e-08 to counter rounding errors.
#' @template treeScorerDots
#' 
#' @return This function returns a tree modified by parsimony ratchet iterations.
#'
#' @references 
#'  \insertRef{Nixon1999}{TreeSearch}
#' 
#' @author Martin R. Smith
#' 
#' @seealso \code{\link{TreeSearch}}
##### @seealso \code{\link{Sectorial}}
#' @seealso Adapted from \code{\link[phangorn]{pratchet}} in the \pkg{phangorn} package.
#' 
#' @keywords  tree 
#' @export
Ratchet <- function (tree, dataset, 
                     InitializeData = PhyDat2Morphy,
                     CleanUpData    = UnloadMorphy,
                     TreeScorer     = MorphyLength,
                     Bootstrapper   = MorphyBootstrap,
                     swappers = list(TBRSwap, SPRSwap, NNISwap),
                     BootstrapSwapper = swappers[[length(swappers)]],
                     returnAll=FALSE, stopAtScore=NULL,
                     ratchIter=100, ratchHits=10, searchIter=2000, searchHits=40,
                     bootstrapIter=ceiling(searchIter/5), bootstrapHits=ceiling(searchHits/5),
                     verbosity=1L, 
                     suboptimal=1e-08, ...) {
  epsilon <- 1e-08
  hits <- 0L
  # initialize tree and data
  if (dim(tree$edge)[1] != 2 * tree$Nnode) stop("tree must be bifurcating; try rooting with ape::root")
  tree <- RenumberTips(tree, names(dataset))
  edgeList <- MatrixToList(tree$edge)
  edgeList <- RenumberEdges(edgeList[[1]], edgeList[[2]])

  initializedData <- InitializeData(dataset)
  on.exit(initializedData <- CleanUpData(initializedData))
  
  bestScore <- if (is.null(attr(tree, 'score'))) {
    TreeScorer(edgeList[[1]], edgeList[[2]], initializedData, ...)
  } else {
    attr(tree, 'score')
  }
  if (verbosity > 0L) cat("\n* Beginning Parsimony Ratchet, with initial score", bestScore)
  if (!is.null(stopAtScore) && bestScore < stopAtScore + epsilon) return(tree)
  if (class(swappers) == 'function') swappers <- list(swappers)
  
  if (returnAll) {
    nullForest <- vector('list', ratchIter)
    forest <- nullForest
    forestScores <- rep(NA, ratchIter)
  }

  iterationsWithBestScore <- 0
  iterationsCompleted <- 0
  for (i in 1:ratchIter) {
    if (verbosity > 1L) cat ("\n* Ratchet iteration", i, "- Generating new tree by bootstrapping dataset. ")
    candidate <- Bootstrapper(edgeList, initializedData, maxIter=bootstrapIter, maxHits=bootstrapHits, 
                                  verbosity=verbosity-2L, EdgeSwapper=BootstrapSwapper, ...)
    candScore <- 1e+08
    if (verbosity > 2L) cat ("\n - Rearranging from new candidate tree:")
    for (EdgeSwapper in swappers) {
      candidate <- EdgeListSearch(candidate, dataset=initializedData, TreeScorer=TreeScorer, 
                                  EdgeSwapper=EdgeSwapper, maxIter=searchIter, maxHits=searchHits,
                                  verbosity=verbosity-2L, ...)                                  
      candScore <- candidate[[3]]
      if (!is.null(stopAtScore) && candScore < stopAtScore + epsilon) break;
    }
    
    if (verbosity > 2L) cat("\n - Rearranged candidate tree scored ", candScore)
    if (returnAll && candScore < (bestScore + suboptimal)) { # Worth saving this tree in forest
      forest[[i]] <- candidate
      forestScores[i] <- candScore
    }
    if ((candScore + epsilon) < bestScore) {
      # New 'best' tree
      edgeList <- candidate
      bestScore <- candScore
      iterationsWithBestScore <- 1L
    } else if (bestScore + epsilon > candScore) { # i.e. best == cand, allowing for floating point error
      iterationsWithBestScore <- iterationsWithBestScore + 1L
      edgeList <- candidate
    }
    if (verbosity > 1L) cat("\n* Best score after", i, "/", ratchIter, "ratchet iterations:",
                            bestScore, "( hit", iterationsWithBestScore, "/", ratchHits, ")")
    if ((!is.null(stopAtScore) && bestScore < stopAtScore + epsilon) 
    || (iterationsWithBestScore >= ratchHits)) {
      iterationsCompleted <- i
      break()
    }
  } # end for
  if (iterationsCompleted == 0) iterationsCompleted <- ratchIter
  if (verbosity > 0L) cat ("\nCompleted parsimony ratchet after", iterationsCompleted, "iterations with score", bestScore, "\n")
   
  if (returnAll) {
    keepers <- !is.na(forestScores) & forestScores < bestScore + suboptimal
    forestScores <- forestScores[keepers]
    forest <- forest[keepers]
    if (verbosity > 1L) cat("\n - Keeping", sum(keepers), "trees from iterations numbered:\n   ", which(keepers))
    if (length(forest) > 1) {
      forest <- lapply(forest, function (phy) {
        x <- tree
        x$edge <- ListToMatrix(phy)
        attr(x, 'score') <- phy[[3]]
        # Return to lapply: 
        x})
      class(forest) <- 'multiPhylo'
      ret <- unique(forest)
      if (verbosity > 1L) cat("\n - Removing duplicates leaves", length(ret), "unique trees")
      uniqueScores <- vapply(ret, attr, double(1), 'score')
    } else if (length(forest) == 1) {
      ret <- tree
      ret$edge <- ListToMatrix(forest[[1]])
      uniqueScores <- forest[[1]][[3]]
    } else {
      stop("\nNo trees!? Is suboptimal set to a sensible (positive) value?")
    }
    if (verbosity > 0L) cat('\nFound', sum(uniqueScores == min(uniqueScores)), 'unique MPTs and', length(ret) - sum(uniqueScores == min(uniqueScores)), 'suboptimal trees.\n')
    # Return:
    ret
  } else {
    tree$edge <- ListToMatrix(edgeList)
    attr(tree, 'score') <- bestScore
    # Return:
    tree  
  }
}

#' @describeIn Ratchet Shortcut for Ratchet search under Profile Parsimony
#' @export
ProfileRatchet <- function (tree, dataset,
                            swappers = list(TBRSwap, SPRSwap, NNISwap),
                            BootstrapSwapper = swappers[[length(swappers)]],
                            returnAll=FALSE, stopAtScore=NULL,
                            ratchIter=100, ratchHits=10, searchIter=2000, searchHits=40,
                            bootstrapIter=searchIter, bootstrapHits=searchHits, verbosity=1L, 
                            suboptimal=1e-08, ...) {
  Ratchet(tree=tree, dataset=dataset,
          InitializeData=ProfileInitMorphy, CleanUpData=ProfileDestroyMorphy,
          TreeScorer=ProfileScoreMorphy, Bootstrapper=ProfileBootstrap,
          swappers=swappers, BootstrapSwapper=BootstrapSwapper,
          returnAll=returnAll, suboptimal=suboptimal, stopAtScore=stopAtScore,
          ratchIter=ratchIter, ratchHits=ratchHits,
          searchIter=searchIter, searchHits=searchHits,
          bootstrapIter=searchIter, bootstrapHits=bootstrapHits, 
          verbosity=verbosity,  ...)
}

#' @describeIn Ratchet Shortcut for Ratchet search using implied weights
#' @template concavityParam
#' @export
IWRatchet <- function (tree, dataset, concavity=4,
                            swappers = list(TBRSwap, SPRSwap, NNISwap),
                            BootstrapSwapper = swappers[[length(swappers)]],
                            returnAll=FALSE, stopAtScore=NULL,
                            ratchIter=100, ratchHits=10, searchIter=2000, searchHits=40,
                            bootstrapIter=searchIter, bootstrapHits=searchHits, verbosity=1L, 
                            suboptimal=1e-08, ...) {
  dataset <- PrepareDataIW(dataset)
  if (verbosity > 1L) cat("\n* Using implied weighting with concavity constant k =", concavity)
  
  Ratchet(tree=tree, dataset=dataset, 
          concavity=concavity, minSteps=attr(dataset, 'min.steps'), 
          InitializeData=IWInitMorphy, CleanUpData=IWDestroyMorphy,
          TreeScorer=IWScoreMorphy, Bootstrapper=IWBootstrap,
          swappers=swappers, BootstrapSwapper=BootstrapSwapper,
          returnAll=returnAll, suboptimal=suboptimal, stopAtScore=stopAtScore,
          ratchIter=ratchIter, ratchHits=ratchHits,
          searchIter=searchIter, searchHits=searchHits,
          bootstrapIter=searchIter, bootstrapHits=bootstrapHits, 
          verbosity=verbosity, ...)
}

#' @describeIn Ratchet returns a list of optimal trees produced by nSearch Ratchet searches
#' @param nSearch Number of Ratchet searches to conduct (for RatchetConsensus)
#' @export
RatchetConsensus <- function (tree, dataset, ratchHits=10, 
                              searchIter=500, searchHits=20, verbosity=0L, 
                              swappers=list(RootedNNISwap), nSearch=10, 
                              stopAtScore=NULL, ...) {
  trees <- lapply(logical(nSearch), function (x) Ratchet(tree, dataset, ratchIter=1, 
                                                         searchIter=searchIter, searchHits=searchHits, verbosity=verbosity, swappers=swappers, 
                                                         stopAtScore=stopAtScore, ...))
  scores <- vapply(trees, function (x) attr(x, 'score'), double(1))
  trees <- unique(trees[scores == min(scores)])
  cat ("Found", length(trees), 'unique trees from', nSearch, 'searches.')
  return (trees)
}


#' @describeIn Ratchet returns a list of optimal trees produced by nSearch 
#'                     Ratchet searches, using implied weighting
#' @export
IWRatchetConsensus <- function (tree, dataset, ratchHits=10, concavity=4,
                              searchIter=500, searchHits=20, verbosity=0L, 
                              swappers=list(RootedNNISwap), nSearch=10, 
                              suboptimal=suboptimal,
                              stopAtScore=NULL, ...) {
  trees <- lapply(logical(nSearch), function (x) IWRatchet(tree, dataset, ratchIter=1, 
                                                         searchIter=searchIter, searchHits=searchHits, verbosity=verbosity, swappers=swappers, 
                                                         stopAtScore=stopAtScore, ...))
  scores <- vapply(trees, function (x) attr(x, 'score'), double(1))
  trees <- unique(trees[scores == min(scores)])
  cat ("Found", length(trees), 'unique trees from', nSearch, 'searches.')
  return (trees)
}
