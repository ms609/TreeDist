#' Rearrange phylogenetic tree
#' @details \code{RearrangeTree} performs one tree rearrangement of a specified type
#' 
#' @template treeParent
#' @template treeChild
#' @param dataset Third argument to pass to \code{TreeScorer}.
#' @template TreeScorerParam
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
#' 
#' @examples
#' data('Lobo')
#' random.tree <- RandomTree(Lobo.phy)
#' RearrangeTree(random.tree, dataset=Lobo.phy, Rearrange=RootedNNI)
#' 
#' @export
RearrangeEdges <- function (parent, child, dataset, TreeScorer, scoreToBeat=TreeScorer(parent, child, dataset),
                            hits=0L, EdgeSwapper, iter='?', verbosity=0L, ...) {
  eps <- 1e-08
  rearrangedEdges <- EdgeSwapper(parent, child) 
  # TODO we probably want to get ALL trees 1 REARRANGE step away
  # One benefit of this is that if NONE of these trees are as good or better, we can give up immediately, 
  # as there's no way out of this local optimum.
  candidateScore <- TreeScorer(rearrangedEdges[[1]], rearrangedEdges[[2]], dataset)
  if (candidateScore > (scoreToBeat + eps)) {
    if (verbosity > 3L) cat("\n    . Iteration", iter, '- Previous score', scoreToBeat, "<", candidateScore)
    return (list(parent, child, scoreToBeat, hits)) # Send back original tree
  } else if (candidateScore + eps > scoreToBeat) { # i.e. scores are equal
    hits <- hits + 1L
    if (verbosity > 2L) cat("\n    - Iteration", iter, "- Best score", scoreToBeat, "hit", hits, "times")
  } else {
    hits <- 1L
    if (verbosity > 1L) cat("\n    * Iteration", iter, "- New best score", candidateScore, "found on", hits, "trees")
  }
  rearrangedEdges[3:4] <- c(candidateScore, hits)
  # Return:
  rearrangedEdges
}


#' Rearrange phylogenetic tree
#' @details \code{MorphyRearrangeTree} performs one tree rearrangement of a specified type
#' 
#' @param tree a rooted bifurcating phylogenetic tree with the desired outgroup, with its labels
#'             in an order that matches the Morphy object, and the attributes
#'             \code{score}, the tree's optimality score, and 
#'             \code{hits}, the number of times the best score has been hit in the calling function.
#' @template morphyObjParam
#' @param Rearrange a rearrangement function that returns a tree: probably one of 
#'     \code{\link{RootedNNI}}, \code{\link{RootedSPR}} or
#'     \code{\link{RootedTBR}}.
#' @param minScore trees longer than \code{min.score}, probably the score of the starting tree,
#'     will be discarded.
#' @param returnSingle returns all trees if \kbd{FALSE} or a randomly selected tree if \kbd{TRUE}.
#' @param iter iteration number of calling function, for reporting to user only.
#' @template clusterParam
#' @template verbosityParam
#' 
#' @return{This function returns the most parsimonious of the trees generated, with attributes \code{hits} and \code{score}
#'  as described for argument \code{tree}, and with tip labels ordered to match morphyObj.}
#' @author Martin R. Smith
#' @seealso
#'   \itemize{
#'     \item \code{\link{RootedNNI}}
#'     \item \code{\link{RootedSPR}}
#'     \item \code{\link{RootedTBR}}
#'   }
#' 
#' @importFrom parallel clusterCall
#' @export
MorphyRearrangeTree <- function (tree, morphyObj, Rearrange, minScore=NULL, returnSingle=TRUE,
                                 iter='?', cluster=NULL, verbosity=0L) {
  if (is.null(attr(tree, 'score'))) bestScore <- 1e+07 else bestScore <- attr(tree, 'score')
  if (is.null(attr(tree, 'hits'))) hits <- 1 else hits <- attr(tree, 'hits')
  if (is.null(cluster)) {
    rearrangedTree <- Rearrange(tree)
    #rearrangedTree <- RenumberTips(Rearrange(tree), tipOrder)
    trees <- list(rearrangedTree)
    minScore <- MorphyTreeLength(rearrangedTree, morphyObj)
    bestTrees <- c(TRUE)
  } else {
    stop("Cluster not implemented.")
    # candidates <- clusterCall(cluster, function(re, tr, k) {ret <- re(tr); attr(ret, 'score') <- Fitch(ret, cl.data, k); ret}, rearrange, tree, concavity)
    # scores <- vapply(candidates, function(x) attr(x, 'ps'), 1)
    # candidates <- lapply(seq_along(cl), function (x) Rearrange(tree)) # TODO don't pick the same tree twice
    # warning("Not tested; likely to fail.")
    # 
    # scores <- parLapply(cluster, seq_along(cluster), function (i) MorphyTreeLength(candidates[[i]], morphyObj[[i]])) # ~3x faster to do this in serial in r233.
    # minScore <- min(scores)
    # bestTrees <- scores == minScore 
    # trees <- candidates[bestTrees]
  }
  if (bestScore < minScore) {
    if (verbosity > 3L) cat("\n    . Iteration", iter, '- Min score', minScore, ">", bestScore)
  } else if (bestScore == minScore) {
    hits <- hits + sum(bestTrees)
    if (verbosity > 2L) cat("\n    - Iteration", iter, "- Best score", minScore, "hit", hits, "times")
  } else {
    hits <- sum(bestTrees)
    if (verbosity > 1L) cat("\n    * Iteration", iter, "- New best score", minScore, "found on", hits, "trees")
  }
  if (length(returnSingle) && returnSingle) trees <- sample(trees, 1)[[1]]
  attr(trees, 'hits') <- hits
  attr(trees, 'score') <- minScore
  trees
}

#' @describeIn MorphyRearrangeTree optimised version that requires parent and child vectors to be extracted from a tree
#' @template treeParent
#' @template treeChild
#' @param inputScore the score of the tree, if known
#' @param hits number of times that this score has been hit
#' @template EdgeSwapperParam
#' @return a rearranged edgeList.
#'
#' @author Martin R. Smith
#' @keywords internal
#' @export
MorphyRearrange <- function (parent, child, morphyObj, inputScore=1e+07, hits=0, 
                             EdgeSwapper, minScore=NULL, returnSingle=TRUE,
                             iter='?', cluster=NULL, verbosity=0L) {
  RearrangeEdges(parent, child, dataset=morphyObj, TreeScorer=MorphyLength, 
                 inputScore=inputScore, hits=hits, EdgeSwapper=EdgeSwapper, 
                 minScore=minScore, returnSingle=returnSingle, iter=iter,
                 cluster=cluster, verbosity=verbosity)
}
