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