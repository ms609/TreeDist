#' Parsimony Ratchet
#'
#' \code{RatchetSearch} uses the parsimony ratchet (Nixon 1999) to search for a more parsimonious tree.
#'
#' @template treeParam 
#' @template datasetParam
#' @param keepAll Set to \kbd{TRUE} to report all MPTs encountered during the search, perhaps to analyze consensus
#' @param maxIt   maximum ratchet iterations to perform;
#' @param maxIter maximum rearrangements to perform on each bootstrap or ratchet iteration;
#' @param maxHits maximum times to hit best score before terminating a tree search within a pratchet iteration;
#' @param k stop when k ratchet iterations have found the same best score;
#' @template stopAtScoreParam
#' @template verbosityParam
#' @param rearrangements (list of) function(s) to use when rearranging trees
#'        e.g. \code{list(TreeSearch::RootedTBR, TreeSearch::RootedNNI)}
#' @param \dots other arguments to pass to subsequent functions.
#' @param nSearch Number of Ratchet searches to conduct (for RatchetConsensus)
#' 
#' @return This function returns a tree modified by parsimony ratchet iteration, retaining the position of the root.
#'
#' @references Nixon, K. C. (1999). \cite{The Parsimony Ratchet, a new method for rapid parsimony analysis.}
#'  Cladistics, 15(4), 407-414. doi:\href{http://dx.doi.org/10.1111/j.1096-0031.1999.tb00277.x}{10.1111/j.1096-0031.1999.tb00277.x}
#'
#' @author Martin R. Smith
#' 
#' Adapted from \code{\link[phangorn]{pratchet}} in the \pkg{phangorn} package, which does not preserve the position of the root.
#' 
#' @seealso \code{\link[phangorn]{pratchet}}
#' @seealso \code{\link{BasicSearch}}
### #' @seealso \code{\link{SectorialSearch}}
#' 
#' @examples{
#' data('inapplicable.datasets')
#' my.phyDat <- inapplicable.phyData[[1]]
#' RatchetSearch(tree=TreeSearch::RandomTree(my.phyDat, root=names(my.phyDat)[1]), 
#'         dataset=my.phyDat, maxIt=1, maxIter=50)
#' }
#' @keywords  tree 
#' @export
RatchetSearch <- function 
(tree, dataset, keepAll=FALSE, maxIt=100, maxIter=5000, maxHits=40, k=10, stopAtScore=NULL,
  verbosity=1L, rearrangements=list(TBRCore, SPRCore, NNICore), ...) {
  if (class(dataset) != 'phyDat') stop("dataset must be of class phyDat, not", class(dataset))
  morphyObj <- LoadMorphy(dataset)
  on.exit(morphyObj <- UnloadMorphy(morphyObj))
  tree <- RenumberTips(tree, names(dataset))
  edge <- tree$edge
  edgeList <- RenumberEdges(edge[, 1], edge[, 2])
  eps <- 1e-08
  bestScore <- attr(tree, "score")
  if (is.null(bestScore)) {
    bestScore <- MorphyLength(edgeList[[1]], edgeList[[2]], morphyObj, ...)
  }
  if (verbosity > 0L) cat("* Initial score:", bestScore)
  if (!is.null(stopAtScore) && bestScore < stopAtScore + eps) return(tree)
  if (keepAll) forest <- vector('list', maxIter)

  if (class(rearrangements) == 'function') rearrangements <- list(rearrangements)
  kmax <- 0 
  for (i in 1:maxIt) {
    if (verbosity > 0L) cat ("\n* Ratchet iteration", i, "- Running NNI on bootstrapped dataset. ")
    candidate <- MorphyBootstrap(edgeList=edgeList, morphyObj=morphyObj, maxIter=maxIter, maxHits=maxHits,
                               verbosity=verbosity-1L, ...)
    
    for (Rearrange in rearrangements) {
      if (verbosity > 0L) cat ("\n - Rearranging new candidate tree...")
      candidate <- DoTreeSearch(candidate, morphyObj, Rearrange=Rearrange, stopAtScore=stopAtScore,
                                verbosity=verbosity-1L, maxIter=maxIter, maxHits=maxHits, ...)
      candScore <- candidate[[3]]
      if (!is.null(stopAtScore) && candScore < stopAtScore + eps) return(candidate)
    }
    if ((candScore + eps) < bestScore) {
      if (keepAll) {
        forest <- vector('list', maxIter)
        forest[[i]] <- candidate
      }
      tree <- candidate
      bestScore <- candScore
      kmax <- 1
    } else {
      if (bestScore + eps > candScore) { # i.e. best == cand, allowing for floating point error
        kmax <- kmax + 1
        candidate$tip.label <- names(dataset)
        tree <- candidate
        if (keepAll) forest[[i]] <- candidate
      }
    }
    if (verbosity > 0L) cat("\n* Best score after", i, "/", maxIt, "ratchet iterations:", 
                           bestScore, "( hit", kmax, "/", k, ")\n")
    if (kmax >= k) break()
  } # for
  if (verbosity > 0L)
    cat ("\n* Completed ratchet search with score", bestScore, "\n")
    
  if (keepAll) {
    forest <- forest[!vapply(forest, is.null, logical(1))]
    class(forest) <- 'multiPhylo'
    ret <- unique(forest)
    cat('Found', length(ret), 'unique MPTs.')
  } else {
    ret <- tree
    ret$tip.label <- names(dataset)
    attr(ret, 'hits') <- NULL
  }
  return (ret)
}

#' @describeIn RatchetSearch returns a list of optimal trees produced by nSearch Ratchet searches
#' @export
RatchetConsensus <- function (tree, dataset, maxIt=5000, maxIter=500, maxHits=20, k=10, verbosity=0L, 
  rearrangements=list(RootedNNI), nSearch=10, ...) {
  trees <- lapply(1:nSearch, function (x) RatchetSearch(tree, dataset, maxIt=maxIt, 
              maxIter=maxIter, maxHits=maxHits, k=1, verbosity=verbosity, rearrangements=rearrangements, ...))
  scores <- vapply(trees, function (x) attr(x, 'score'), double(1))
  trees <- unique(trees[scores == min(scores)])
  cat ("Found", length(trees), 'unique trees from ', nSearch, 'searches.')
  return (trees)
}

#' Bootstrap tree search with inapplicable data
#' 
#' @template edgeListParam
#' @template morphyObjParam
#' @param maxIter maximum number of iterations to perform in tree search
#' @param maxHits maximum number of hits to accomplish in tree search
#' @template verbosityParam
#' @param \dots further parameters to send to \code{DoTreeSearch}
#'
#' @return A tree that is optimal under a random sampling of the original characters
#' @export
MorphyBootstrap <- function (edgeList, morphyObj, maxIter, maxHits, verbosity=1L, ...) {
## Simplified version of phangorn::bootstrap.phyDat, with bs=1 and multicore=FALSE
  startWeights <- MorphyWeights(morphyObj)[1, ]
  eachChar <- seq_along(startWeights)
  v <- rep(eachChar, startWeights)
  BS <- tabulate(sample(v, replace=TRUE), length(startWeights))
  vapply(eachChar, function (i) 
         mpl_set_charac_weight(i, BS[i], morphyObj), integer(1))
  mpl_apply_tipdata(morphyObj)
  # TODO Optimisation exploration: is NNICore better than TBRCore?
  res <- DoTreeSearch(edgeList, morphyObj, Rearrange=NNICore, maxIter=maxIter, maxHits=maxHits, verbosity=verbosity-1L, ...)
  vapply(eachChar, function (i) 
         mpl_set_charac_weight(i, startWeights[i], morphyObj), integer(1))
  mpl_apply_tipdata(morphyObj)
  res[1:2]
}

#' DoTreeSearch
#'
#' Performs a tree search
#' 
#' Does the hard work of searching for a most parsimonious tree, given the
#' parent and child vectors of a tree arranged in preorder (perhaps with 
#' \code{\link[TreeSearch]{RenumberEdges}}).
#' End-users are expected to access this function through its wrapper, TreeSearch
#' It is also called directly by RatchetSearch and Sectorial functions
#'
#' @template edgeListParam
#' @template morphyObjParam
#' @param Rearrange Function to use to rearrange trees; example: 
#'                  \code{TreeSearch::\link[TreeSearch]{RootedTBR}}.
#' @param maxIter maximum iterations to conduct.
#' @param maxHits stop search after this many hits.
#' @template stopAtScoreParam
#' @param forestSize how many trees to hold.
#' @template clusterParam
#' @template verbosityParam
#' @param \dots additional variables to pass to \code{\link{MorphyRearrangeTree}}.
#'
#' @author Martin R. Smith
#' 
#' @keywords internal
#' @export

DoTreeSearch <- function (edgeList, morphyObj, Rearrange, maxIter=100, maxHits=20, 
                          stopAtScore=NULL, forestSize=1L, 
                          cluster=NULL, verbosity=1L, ...) {
  eps <- 1e-07                        
  if (!is.null(forestSize) && length(forestSize)) {
    if (forestSize > 1L) {
      stop("TODO: Forests not supported")
      forest <- empty.forest <- vector('list', forestSize)
      forest[[1]] <- tree
    } else {
      forestSize <- 1L
    }
  }
  if (length(edgeList) < 3) {
    bestScore <- MorphyLength(edgeList[[1]], edgeList[[2]], morphyObj)
  } else {
    bestScore <- edgeList[[3]]
  }
  hits <- if (length(edgeList) < 4) 0 else edgeList[[4]]
  if (verbosity > 0L) cat("  - Initial score:", bestScore)
  if (!is.null(stopAtScore) && bestScore < stopAtScore + eps) return(edgeList)
  returnSingle <- !(forestSize > 1L)
  
  for (iter in 1:maxIter) {
    candidateLists <- MorphyRearrange(edgeList[[1]], edgeList[[2]], morphyObj, 
                             inputScore=bestScore, hits=hits, RearrangeEdges=Rearrange,
                             minScore=bestScore, returnSingle=returnSingle, iter=iter, 
                             cluster=cluster, verbosity=verbosity, ...)
    scoreThisIteration <- candidateLists[[3]]
    if (forestSize > 1L) {
      stop("TODO re-code this")
      if (scoreThisIteration == bestScore) {
        forest[(hits-length(edgeLists)+1L):hits] <- trees ## TODO Check that length still hojlds
        edgeList <- sample(forest[1:hits], 1)[[1]]
        attr(edgeList, 'score') <- scoreThisIteration
        attr(edgeList, 'hits') <- hits
      } else if (scoreThisIteration < bestScore) {
        bestScore <- scoreThisIteration
        forest <- empty.forest
        forest[1:hits] <- candidateLists
        edgeList <- sample(candidateLists , 1)[[1]]
        attr(edgeList, 'score') <- scoreThisIteration
      }
    } else {
      if (scoreThisIteration <= bestScore) {
        bestScore <- scoreThisIteration
        edgeList <- candidateLists
        if (!is.null(stopAtScore) && bestScore < stopAtScore + eps) return(edgeList)
      }
    }
    if (hits >= maxHits) break
  }
  if (verbosity > 0L) cat("\n  - Final score", attr(edgeList, 'score'), "found", attr(edgeList, 'hits'), "times after", iter, "rearrangements\n")  
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
#' @param tree a fully-resolved starting tree in \code{\link{phylo}} format, with the desired outgroup; 
#'        edge lengths are not supported and will be deleted.
#' @template datasetParam
#' @param Rearrange Function used to rearrange trees; default: \code{\link[TreeSearch]{RootedTBR}}.
#' @param maxIter the maximum number of iterations to perform before abandoning the search.
#' @param maxHits the maximum times to hit the best score before abandoning the search.
#' @param forestSize the maximum number of trees to return - useful in concert with \code{\link{consensus}}.
#' @template nCoresParam
#' @template verbosityParam
#' @param \dots other arguments to pass to subsequent functions.
#' 
#' @return{
#' This function returns a tree, with an attribute \code{score} conveying its parsimony score.
#' Note that the parsimony score will be inherited from the tree's attributes, which is only valid if it 
#' was generated using the same \code{data} that is passed here.
#' }
#' @author Martin R. Smith
#'
#' @seealso
#' \itemize{
#' \item \code{\link{InapplicableFitch}}, calculates parsimony score, supports inapplicable tokens;
#' \item \code{\link[TreeSearch]{RootedNNI}}, conducts tree rearrangements;
### #' \item \code{\link{SectorialSearch}}, alternative heuristic, useful for larger trees;
#' \item \code{\link{RatchetSearch}}, alternative heuristic, useful to escape local optima.
#' }
#'
#' @examples
#' data('inapplicable.datasets')
#' my.phyDat <- inapplicable.phyData[[1]]
#' njtree <- TreeSearch::NJTree(my.phyDat)
#'
#' \dontrun{
#' TreeSearch(njtree, my.phyDat, maxIter=20, Rearrange=TreeSearch::NNI)
#' TreeSearch(njtree, my.phyDat, maxIter=20, Rearrange=TreeSearch::RootedSPR)
#' }
#' 
#' @keywords  tree 
#' 
#' @export
BasicSearch <- function 
(tree, dataset, Rearrange=TBRCore, maxIter=100, maxHits=20, forestSize=1, 
 nCores=1L, verbosity=1, ...) {
  # Initialize morphy object
  if (class(dataset) != 'phyDat') stop ("dataset must be of class phyDat, not ", class(dataset))
  if (dim(tree$edge)[1] != 2 * tree$Nnode) stop("tree must be bifurcating; try rooting with ape::root")
  tree <- RenumberTips(tree, names(dataset))
  edge <- tree$edge
  edgeList <- RenumberEdges(edge[, 1], edge[, 2])
  if (nCores > 1L) {
    stop("Clusters are not yet supported (#23).")
    ### cluster <- snow::makeCluster(nCores)
    ### on.exit(snow::stopCluster(cluster), add=TRUE)
    ### snow::clusterEvalQ(cluster, {library(inapplicable); NULL})
    ### morphyObj <- lapply(seq_len(nCores), function(xx) LoadMorphy(dataset))
    ### on.exit(morphyObj <- vapply(morphyObj, UnloadMorphy, integer(1)), add=TRUE)
    ### snow::clusterExport(cluster, c('dataset'))
  } else {
    morphyObj <- LoadMorphy(dataset)
    cluster <- NULL
    on.exit(morphyObj <- UnloadMorphy(morphyObj))
  }
  ret <- DoTreeSearch(edgeList, morphyObj, Rearrange=Rearrange,
                      maxIter=maxIter, maxHits=maxHits, forestSize=forestSize, cluster=cluster, 
                      verbosity, ...)
  return (ret)
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
###   #' @seealso \code{\link{BasicSearch}}
###   #' @seealso \code{\link{RatchetSearch}}
###   #' 
###   #' @examples
###   #' require('ape')
###   #' data('SigSut')
###   #' outgroup <- c('Lingula', 'Mickwitzia', 'Neocrania')
###   #' njtree <- ape::root(nj(dist.hamming(SigSut.phy)), outgroup, resolve.root=TRUE)
###   #' njtree$edge.length <- NULL; njtree<-ape::root(njtree, outgroup, resolve.root=TRUE)
###   #' InapplicableSectorial(njtree, SigSut.phy, maxIt=1, maxIter=50, largest.sector=7)
###   #' \dontrun{SectorialSearch(njtree, SigSut.phy) # Will be time-consuming }
###   #' 
###   #' 
###   #' @keywords  tree 
###   #' @export
###   SectorialSearch <- function
###   (tree, dataset, SectorialRearrangements=TreeSearch::NNI, maxIter=2000,
###    subsequentRearrangements = list(TreeSearch::RootedNNI, TreeSearch::RootedTBR, 
###     TreeSearch::RootedSPR, TreeSearch::RootedNNI), verbosity=3, ...) {
###     if (class(dataset) != 'phyDat') stop("dataset must be of class phyDat, not", class(dataset))
###     morphyObj <- LoadMorphy(dataset)
###     on.exit(morphyObj <- UnloadMorphy(morphyObj))
###     tree <- TreeSearch::RenumberTips(TreeSearch::Renumber(tree), names(dataset))
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
###     tree <- TreeSearch::RenumberTips(TreeSearch::Renumber(tree), names(dataset))
###     if (length(bestScore) == 0) bestScore <- InapplicableFitch(tree, dataset, ...)[[1]]
###     sect <- MorphySectorial(tree, morphyObj, verbosity=verbosity-1, maxIt=30, 
###       maxIter=maxIter, maxHits=15, smallest.sector=6, 
###       largest.sector=length(tree$edge[,2L])*0.25, rearrangements=rearrangements)
###     sect <- BasicSearch(sect, dataset, Rearrange=subsequentRearrangements, maxIter=maxIter, maxHits=30, cluster=cluster, verbosity=verbosity)
###     if (attr(sect, 'score') <= bestScore) {
###       return (sect)
###     } else return (tree)
###   }
###   