#' Sectorial Search
#'
#' \code{SectorialSearch} performs a sectorial search on a tree, preserving the position of the root.
#'
#' \code{DoSectorial} performs a sectorial search on the tree specified. A sectorial search 
#' detaches a random part of the tree, performs rearrangments on this subtree, then reattaches it 
#' to the main tree (Goloboff, 1999).
#' The improvement to local \var{pscore} hopefully (but not necessarily) improves the overall \var{pscore}.
#' As such, the output of \code{DoSectorial} should be treated by further \acronym{TBR (/SPR/NNI)}
#' rearrangements and only retained if the ultimate parsimony score is better than 
#' that of the original tree.
#' 
#' \code{SectorialSearch} is a basic recipe that runs \code{DoSectorial} followed by a few rounds
#' of tree rearrangement, returning a tree whose \var{pscore} is no worse than that of \code{start.tree}.
#' 
#' @param tree a rooted, resolved tree in \code{\link{phylo}} format from which to start the search;
#' @template datasetParam
#' @param outgroup a vector listing the taxa that form the outgroup;
#' @template concavityParam
#' @param maxIter maximum number of rearrangements to perform on each search iteration - provide a list or vector 
#'                with each entry corresponding to an entry in Rearrangements, or a single entry to be used for all
#' @param maxHits maximum number of hits to accomplish on each search iteration - provide a list or vector 
#'                with each entry corresponding to an entry in Rearrangements, or a single entry to be used for all
#' @param cluster a cluster prepared using \code{\link{PrepareCluster}}; may speed up search on multicore machines;
#' @param k stop when \code{k} searches have improved their sectorial score;
#' @template verbosityParam
#' @param smallestSector sectors with fewer than \code{smallestSector} taxa will not be selected; \kbd{4} is the smallest sensible value;
#' @param largestSector sectors with more than \code{largestSector} taxa will not be selected;
#' @param sectRearrangements A list of functions to be applied when rearranging trees in sectorial 
#'                           search. Note: these must retain the position of the root of the tree
#' @param searchRearrangements A list of functions to be applied, in turn, during tree search.
#'                             Default: \code{list(\link{RootedTBR}, \link{RootedNNI})}
#' @template treeScorerDots
#' 
#' @return a rooted tree of class \code{phylo}.
#' 
#' @references Goloboff, P. (1999). \cite{Analyzing large data sets in reasonable times: solutions for composite optima.}
#' Cladistics, 15(4), 415-428. doi:\href{http://dx.doi.org/10.1006/clad.1999.0122}{10.1006/clad.1999.0122}
#' 
#' @author Martin R. Smith
#' 
#' @seealso \code{\link{TreeSearch}}
#' @seealso \code{\link{Ratchet}}
#' 
#' @examples
#' data('Lobo')
#' startTree <- RandomTree(Lobo.phy, 'Cricocosmia') # Position of root will be fixed
#' firstEstimate <- SectorialSearch(startTree, Lobo.phy, maxIter=50, largestSector=8, verbosity=5)
#' \dontrun{SectorialSearch(firstEstimate, Lobo.phy, sectRearrangements=list(RootedNNI, RootedTBR, RootedNNI)) # Will be time-consuming}
#' 
#' @keywords  tree 
#' @export
SectorialSearch <- function (tree, dataset, TreeScorer = FitchScore, sectRearrangements=list(RootedNNI), 
                             searchRearrangements=list(RootedNNI, RootedTBR, RootedNNI), 
                             maxHits=c(30, 40, 60), maxIter=2000, cluster=NULL, verbosity=3, ...) {
  best.score <- attr(tree, 'score')
  if (is.null(treeOrder <- attr(tree, 'order')) || treeOrder != 'preorder') tree <- Preorder(tree)
 
  tree <- RenumberTips(tree, names(dataset))
  if (length(best.score) == 0) best.score <- TreeScorer(tree, dataset, ...)
  sect <- DoSectorial(tree, dataset, cluster=cluster,
    verbosity=verbosity-1, maxit=30, maxIter=max(maxIter), maxHits=15, smallestSector=6, 
    largestSector=dim(tree$edge)[1]*0.25, Rearrangements=sectRearrangements)
  for (i in seq_along(Rearrangements)) {
    iters <- if (length(maxIter) <= i) maxIter[[i]] else min(maxIter)
    hits  <- if (length(maxHits) <= i) maxHits[[i]] else min(maxHits)
    sect <- TreeSearch(sect, dataset, TreeScorer, Rearrangements[[i]], maxIter=iters,
                         maxHits=hits, cluster=cluster, verbosity=verbosity-1)
  }
  if (attr(sect, 'score') <= best.score) {
    return (sect)
  } else return (tree)
}

#' Sector Data
#' Check that chosen sector contains parsimony-informative data
#'
#' The function simply checks that some characters have more than one state.
#' It's crude, but the cost of a false positive is low.
#'
#' @param dataset the dataset to subsample
#' @param tips character vector listing tips that exist in the sector 
#'
#' @keywords internal
#' @export
SectorHasData <- function (dataset, tips) {
  if (class(dataset) =='phyDat') {
    levs <- attr(dataset, 'levels')
    contrast <- attr(dataset, 'contrast')
    index <- as.integer(contrast %*% 2L ^ (seq_along(attr(dataset, 'levels')) - 1))
    characters <- vapply(dataset, function (X) index[X], integer(attr(dataset, 'nr')))
  } else if (names(dataset)) {
    characters <- vapply(dataset, c, integer(length(dataset[[1]])))
  } else if(rownames(dataset)) {
    characters <- t(dataset)
  } else characters <- dataset
  tokens <- apply(characters[, tips], 2, function (x) length(unique(x)))
  return (any(tokens > 1))
}

#' Sectorial Search with inapplicable data
#'
#' \code{DoSectorial} is called by SectorialSearch
#'
#' @template preorderTreeParam
#' @param dataset a dataset in the format expected by \code{TreeScorer}
#' @param TreeScorer a function that will score a tree topology
#' @param maxSectIter maximum number of sectorial iterations to perform
#' @param maxIter maximum number of iterations to perform in tree rearrangement functions
#' @param maxImprovements maximum number of times to find an optimal score before ending sectorial search
#' @param Rearrangements a list of tree rearrangement functions that retain the root of the tree
#'                       (e.g. \code{list(RootedSPR)})
#' @param smallestSector integer giving size of smallest sector to rearrange
#' @param largestSector integer giving size of largest sector to rearrange (rounded down if non-integral)
#' @template verbosityParam
#' @template treeScorerDots
#'  
#' @return a tree of class \code{phylo} with a \code{TreeScorer} score as good or better than that of \code{tree}
#'   
#' @author Martin R. Smith
#' @export
DoSectorial <- function (tree, dataset, TreeScorer = FitchScore, maxSectIter=100, 
                         maxIter=500, maxImprovements=5, smallestSector=4, largestSector=1e+06, 
                         Rearrangements=list(RootedNNI), verbosity=0, ...) {
  tipOrder <- names(dataset)
  if (is.null(tipOrder)) tipOrder <- colnames(dataset)
  if (is.null(tipOrder)) tipOrder <- rownames(dataset)
  if (is.null(tipOrder)) stop("Could not find tips in dataset")
  tree <- RenumberTips(tree, tipOrder)
  if (verbosity >= 0) cat('DoSectorial search: optimizing sectors of', smallestSector, 'to', floor(largestSector), 'tips')
  nEdge <- dim(tree$edge)[1]
  nTip <- (nEdge / 2) + 1
  nonRootNodes <- (nTip + 2):(nEdge + 1)
  
  eps <- 1e-08
  improvements <- 1
  for (i in seq_len(maxSectIter)) {
    nodeLengths <- CladeSizes(tree, nonRootNodes)
    candidateNodes <- nonRootNodes[nodeLengths >= smallestSector & nodeLengths <= largestSector]
    if (verbosity >= 0) cat ("\n - Iteration", i, "- attempting sectorial search on node ")
    repeat {
      sector <- sample(candidateNodes, 1)
      candidate <- Subtree(tree, sector)
      crownTips <- candidate$tip.label
      sectorSize <- length(crownTips)
      cat(sector, 'size', sectorSize, '...')
      
      if (SectorHasData(dataset, crownTips)) break else cat('unsuitable (no dataset); trying')
      
      candidateNodes <- candidateNodes[-match(sector, candidateNodes)]
      if (length(candidateNodes == 0)) stop('No selectable sectors contain parsimony information! Either "largestSector" is close to "smallestSector" or your dataset is short of parsimony information.')
    }
    if (verbosity >= 0) cat(' Sector OK.')
    ###
    #####crown <- root(AddTip(crown, 0, 'SECTOR_ROOT'), length(crown$tip.label) + 1, resolve.root=TRUE) ## TODO use Root or add ape::root to includeFrom in NAMESPACE
    initialScore <- TreeScorer(candidate, dataset, ...)
    attr(candidate, 'score') <- initialScore
    
    if (verbosity >= 0) cat("\n - Rearranging sector", sector)
    for (Rearrange in Rearrangements) {
      candidate <- TreeSearch(candidate, dataset, TreeScorer, Rearrange,
                                verbosity=verbosity-1, maxIter=maxIter, ...) 
    }
    candidateScore <- attr(candidate, 'score')
    
    if((candidateScore + eps) < initialScore) {
      improvements <- improvements + 1
      
      subtree.labels <- crownTips
      subtree.nTips <- sectorSize
      subtree.edge <- candidate$edge
      subtree.parent <- subtree.edge[, 1]
      subtree.child  <- subtree.edge[, 2]
           
      isTip <- subtree.child <= subtree.nTips
      subtree.child[isTip] <- match(crownTips[subtree.child[isTip]], tipOrder)
      nodeAdjust <- sector - (subtree.nTips + 1)

      subtree.child[!isTip] <- subtree.child[!isTip] + nodeAdjust
      edges <- which(tree$edge[, 2] == sector) + seq_along(subtree.parent)
      tree$edge[edges, 1] <- subtree.parent + nodeAdjust
      tree$edge[edges, 2] <- subtree.child
      
      if (verbosity > 0) cat(' : improved local pscore, updated tree')
    } else if (verbosity > 0) cat (' : no improvement to local pscore')
    if (improvements == maxImprovements) break()
  } # for
  if (verbosity >= 0)
    cat ("\nCompleted sectorial rearrangements.\n")
  attr(tree, 'score') <- NULL
  attr(tree, 'hits') <- NULL
  tree
}  # DoSectorial

