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
  ret <- DoTreeSearch(tree, dataset, method, maxiter, maxhits, forest.size, cluster, 
                      verbosity, ...)
  return (ret)
}

#' DoTreeSearch
#'
#' Performs a tree search
#' 
#' Does the hard work of searching for a most parsimonious tree.
#' End-users are expected to access this function through its wrapper, TreeSearch
#' It is also called directly by Ratchet and Sectorial functions
#'
#' @template labelledTreeParam
#'
#' @author Martin R. Smith
#' 
#' @keywords internal
#' @export
DoTreeSearch <- function (tree, data, ParsimonyScorer = phangorn:::fitch,  method = 'NNI', 
                        maxiter = 100, maxhits = 20, forest.size = 1,
                        cluster = NULL, track = 1, ...) {
  tree$edge.length <- NULL # Edge lengths are not supported
  attr(tree, 'hits') <- 1
  if (exists("forest.size") && length(forest.size) && forest.size > 1) {
    forest <- empty.forest <- vector('list', forest.size)
    forest[[1]] <- tree
  } else {
    forest.size <- 1 
  }
  if (is.null(attr(tree, 'score'))) attr(tree, 'score') <- ParsimonyScorer(tree, data)
  best.score <- attr(tree, 'score')
  if (track > 0) cat("\n  - Performing", method, "search.  Initial score:", best.score)
  Rearrange <- switch(method, 'TBR' = TBR, 'SPR' = SPR, 'NNI' = NNI)
  return.single <- !(forest.size > 1)
  
  for (iter in 1:maxiter) {
    trees <- RearrangeTree(tree, data, Rearrange, ParsimonyScorer, min.score=best.score,
                           return.single=return.single, iter=iter, cluster=cluster, track=track)
    iter.score <- attr(trees, 'score')
    if (length(forest.size) && forest.size > 1) {
      hits <- attr(trees, 'hits')
      if (iter.score == best.score) {
        forest[(hits-length(trees)+1L):hits] <- trees
        tree <- sample(forest[1:hits], 1)[[1]]
        attr(tree, 'score') <- iter.score
        attr(tree, 'hits') <- hits
      } else if (iter.score < best.score) {
        best.score <- iter.score
        forest <- empty.forest
        forest[1:hits] <- trees
        tree <- sample(trees, 1)[[1]]
        attr(tree, 'score') <- iter.score
        attr(tree, 'hits') <- hits
      }      
    } else {
      if (iter.score <= best.score) {
        best.score <- iter.score
        tree <- trees
      }
    }
    if (attr(trees, 'hits') >= maxhits) break
  }
  if (track > 0) cat("\n  - Final score", attr(tree, 'score'), "found", attr(tree, 'hits'), "times after", iter, method, "iterations\n")  
  if (forest.size > 1) {
    if (hits < forest.size) forest <- forest[-((hits+1):forest.size)]
    attr(forest, 'hits') <- hits
    attr(forest, 'score') <- best.score
    return(unique(forest))
  } else {
    return(tree)
  }
}