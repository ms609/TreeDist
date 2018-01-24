#' Rearrange phylogenetic tree
#' @details \code{RearrangeTree} performs one tree rearrangement of a specified type
#' 
#' @template treeParent
#' @template treeChild
#' @param dataset Third argument to pass to \code{TreeScorer}.
#' @template treeScorerParam
#' @param scoreToBeat Double giving score of input tree.
#' @param hits Integer giving number of times the input tree has already been hit.
#' @template EdgeSwapperParam
#' @param  minScore trees longer than \code{minScore}, probably the score of the best previously known tree,
#'     will be discarded;
##' @param returnSingle returns all trees if \kbd{FALSE} or a randomly selected tree if \kbd{TRUE}.
#' @param iter iteration number of calling function, for reporting to user only.
##' @template clusterParam
#' @template verbosityParam
#' @template treeScorerDots
#'
#' @author Martin R. Smith
#'
#' @template returnEdgeList
## Score will be correct.  Hits will be the number of times it's been hit in this function.
#' 
#' @examples
#' data('Lobo')
#' random.tree <- RandomTree(Lobo.phy)
#' RearrangeTree(random.tree, dataset=Lobo.phy, EdgeSwapper=RootedNNISwap)
#' 
#' @export
RearrangeEdges <- function (parent, child, dataset, TreeScorer, scoreToBeat=TreeScorer(parent, child, dataset),
                            EdgeSwapper, iter='?', hits=0L, verbosity=0L, ...) {
  eps <- 1e-08
  rearrangedEdges <- EdgeSwapper(parent, child)
  # TODO we probably want to get ALL trees 1 REARRANGE step away
  # One benefit of this is that if NONE of these trees are as good or better, we can give up immediately, 
  # as there's no way out of this local optimum.
  candidateScore <- TreeScorer(rearrangedEdges[[1]], rearrangedEdges[[2]], dataset)
  if (candidateScore > (scoreToBeat + eps)) {
    if (verbosity > 3L) cat("\n    . Iteration", iter, '- Rearranged tree score', candidateScore, "> target", scoreToBeat)
  } else if (candidateScore + eps > scoreToBeat) { # i.e. scores are equal
    hits <- hits + 1L
    if (verbosity > 2L) cat("\n    - Iteration", iter, "- Best score", scoreToBeat, "hit", hits, "times")
  } else {
    hits <- 1L
    if (verbosity > 1L) cat("\n    * Iteration", iter, "- New best score", candidateScore, "found on", hits, "trees")
  }
  # TODO when we search multiple trees at once in this function, update hits to the number of best trees.
  rearrangedEdges[3:4] <- c(candidateScore, hits)
  # Return:
  rearrangedEdges
}
