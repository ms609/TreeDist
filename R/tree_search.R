Suboptimality <- function (trees, proportional = FALSE) {
  scores <- vapply(trees, attr, double(1), 'score')
  if (proportional) {
    return ((scores - min(scores)) / min(scores))
  } else {
    return(scores - min(scores))
  }
}

ProfileScore <- function (tree, data) {
    # Data
  if (class(data) == 'phyDat') data <- PrepareDataFitch(data)
  if (class(data) != 'fitchDat') stop('Invalid data type; try ProfileScore(tree, data <- 
                                      PrepareDataFitch(valid.phyDat.object)).')
  at <- attributes(data)
  n.char  <- at$nr # strictly, transformation series patterns; these'll be upweighted later
  weight <- at$weight
  steps <- Fitch(tree, data, at)
  info <- at$info.amounts
  # Return a negative rather than positive value because algorithms assume that 
  # smaller numbers are better
  return(-sum(info[max(0, (steps - 1)) * n.char + seq_len(n.char)] * weight))
}

SuccessiveWeights <- function(tree, data) {
  # Data
  if (class(data) == 'phyDat') data <- PrepareDataFitch(data)
  if (class(data) != 'fitchDat') {
    stop('Invalid data type; prepare data with PhyDat() or PrepareDataFitch().')
  }
  at <- attributes(data)
  weight <- at$weight
  sa.weights <- at$sa.weights
  if (is.null(sa.weights)) sa.weights <- rep(1, length(weight))
  steps <- Fitch(tree, data, at)
  return(sum(steps * sa.weights * weight))
}

SuccessiveApproximations <- function (tree, data, outgroup = NULL, k = 3, max.succiter = 20,
                                      pratchhits = 100, searchhits = 50, searchiter = 500,
                                      pratchiter = 5000, track = 0, suboptimal = 0.1) {
  
  if (k < 1) stop ('k should be at least 1, see Farris 1969 p.379')
  attr(data, 'sa.weights') <- rep(1, length(attr(data, 'weight')))
  collect.suboptimal <- suboptimal > 0
  
  max.node <- max(tree$edge[, 1])
  n.tip <- length(tree$tip.label)
  n.node <- max.node - n.tip
  bests <- vector('list', max.succiter + 1)
  bests.consensus <- vector('list', max.succiter + 1)
  best <- bests[[1]] <- bests.consensus[[1]] <- Root(tree, outgroup)
  for (i in seq_len(max.succiter) + 1) {
    if (track > 0) cat('\nSuccessive Approximations Iteration', i - 1)
    attr(best, 'score') <- NULL
    if (suboptimal > 0) {
      suboptimal.search <- suboptimal * sum(attr(data, 'sa.weights') * attr(data, 'weight'))
    }
    trees <- Ratchet(best, data, ParsimonyScorer = SuccessiveWeights, all = collect.suboptimal, 
                           suboptimal=suboptimal.search,    rearrangements='NNI',
                           pratchhits=pratchhits, searchhits=searchhits, searchiter=searchiter, 
                           pratchiter=pratchiter, outgroup = outgroup, track=track - 1)
    trees <- unique(trees)
    bests[[i]] <- trees
    suboptimality <- Suboptimality(trees)
    bests.consensus[[i]] <- consensus(trees[suboptimality == 0])
    if (all.equal(bests.consensus[[i]], bests.consensus[[i - 1]])) return(bests[2:i])
    best <- trees[suboptimality == 0][[1]]
    l.i <- Fitch(best, data)
    p.i <- l.i / (n.node - 1)
    w.i <- ((p.i)^-k) - 1
    attr(data, 'sa.weights') <- w.i
  }
  cat('Stability not reached.')
  return(bests)
}

#' Search for most parsimonious trees
#'
#' Run standard search algorithms (\acronym{NNI}, \acronym{SPR} or \acronym{TBR}) 
#' to search for a more parsimonious tree.
#'  
#' @param tree a fully-resolved starting tree in \code{\link{phylo}} format, with the desired outgroup; edge lengths are not supported and will be deleted;
#' @template datasetParam
#' @param outgroup a vector listing the taxa in the outgroup;
#' @param concavity concavity constant for implied weighting (not currently implemented!); 
#' @param method rearrangements to perform; one of \kbd{NNI}, \kbd{SPR}, or \kbd{TBR};
#' @param maxiter the maximum number of iterations to perform before abandoning the search;
#' @param maxhits the maximum times to hit the best pscore before abandoning the search;
#' @param forest.size the maximum number of trees to return - useful in concert with \code{\link{consensus}};
#' @param cluster a cluster prepared using \code{\link{PrepareCluster}}; may speed up search on multicore machines;
#' @param verbosity higher values provide more verbose user feedback in stdout;
#' @param \dots other arguments to pass to subsequent functions.
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
#' \item \code{\link{InapplicableFitch}}, calculates parsimony score, supports inapplicable tokens;
#' \item \code{\link{RootedNNI}}, conducts tree rearrangements;
#' \item \code{\link{SectorialSearch}}, alternative heuristic, useful for larger trees;
#' \item \code{\link{Ratchet}}, alternative heuristic, useful to escape local optima.
#' }
#'
#' @examples
#' data('SigSut')
#' outgroup <- c('Lingula', 'Mickwitzia', 'Neocrania')
#' njtree <- Root(nj(dist.hamming(SigSut.phy)), outgroup, resolve.root=TRUE)
#' njtree$edge.length <- NULL; njtree<-SetOutgroup(njtree, outgroup)
#'
#' \dontrun{
#' TreeSearch(njtree, SigSut.phy, outgroup, maxiter=20, method='NNI')
#' TreeSearch(njtree, SigSut.phy, outgroup, maxiter=20, method='SPR')
#' TreeSearch(njtree, SigSut.phy, outgroup, maxiter=20, method='TBR')}
#' 
#' @keywords  tree 
#' 
#' @export
TreeSearch <- function 
(tree, dataset, method='NNI', maxiter=100, maxhits=20, forest.size=1, cluster=NULL, 
 verbosity=1, ...) {
  # Initialize morphy object
  if (class(dataset) != 'phyDat') stop ("dataset must be of class phyDat, not ", class(dataset))
  tree <- RenumberTips(tree, names(dataset))
  morphyObj <- LoadMorphy(dataset)
  ret <- DoTreeSearch(tree, morphyObj, method, maxiter, maxhits, forest.size, cluster, 
                      verbosity, ...)
  morphyObj <- UnloadMorphy(morphyObj)
  return (ret)
}

