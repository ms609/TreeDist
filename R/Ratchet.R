#' Parsimony Ratchet
#'
#' \code{Ratchet} uses the parsimony ratchet (Nixon 1999) to search for a more parsimonious tree.
#'
#' @template treeParam 
#' @param dataset a dataset in the format required by TreeScorer
#' @param TreeScorer Function to score trees, to receive arguments \code{tree} and \code{dataset}.
#' @param returnAll Set to \code{TRUE} to report all MPTs encountered during the search, perhaps to analyze consensus
#' @param rooted whether to retain the position of the root in tree search (TRUE by default)
#' @param ratchIter stop when this many ratchet iterations have been performed
#' @param ratchHits stop when this many ratchet iterations have found the same best score
#' @param searchIter maximum rearrangements to perform on each bootstrap or ratchet iteration;
#' @param searchHits maximum times to hit best score before terminating a tree search within a ratchet iteration;
#' @param bootstrapIter maximum rearrangements to perform on each bootstrap  iteration (default: \code{maxIter})
#' @param bootstrapHits maximum times to hit best score on each bootstrap  iteration (default: \code{maxHits})
#' @template verbosityParam
#' @param rearrangements method(s) to use when rearranging trees. Provide:
#'        the character string "TBR" to conduct TBR, SPR and NNI rearrangements, or "TBR only"
#'        to just perform TBR rearrangements (retaining the position of the root if \code{outgroup = TRUE})
#'        OR: a list of functions to use, one at a time, as the \code{Rearrange} parameter
#'            in successive calls to TreeSearch
#' @param suboptimal retain trees that are suboptimal by this score. Defaults to 1e-08 to counter rounding errors.
#'
#' @template treeScorerDots
#' 
#' @return This function returns a tree modified by parsimony ratchet iteration, retaining the position of the root.
#'
#' @references Nixon, K. C. (1999). \cite{The Parsimony Ratchet, a new method for rapid parsimony analysis.} Cladistics, 15(4), 407-414. doi:\href{http://dx.doi.org/10.1111/j.1096-0031.1999.tb00277.x}{10.1111/j.1096-0031.1999.tb00277.x}
#'
#' @author Martin R. Smith
#' 
#' Adapted from \code{\link[phangorn]{pratchet}} in the \pkg{phangorn} package, which does not preserve the position of the root.
#' 
#' @seealso \code{\link{Ratchet}}
#' @seealso \code{\link{TreeSearch}}
#' @seealso \code{\link{Sectorial}}
#' 
#' @examples Ratchet(RandomTree(Lobo.phy), Lobo.phy, rooted=TRUE,
#'                   ratchHits=2, ratchIter=5, searchIter=200, searchHits=5)
#' 
#' @keywords  tree 
#' @export
## TODO use Rooted NNI / SPR / TBR 
Ratchet <- function (tree, dataset, TreeScorer=FitchScore, returnAll=FALSE, rooted=TRUE, 
                      ratchIter=100, ratchHits=10, searchIter=2000, searchHits=40,
                      bootstrapIter=searchIter, bootstrapHits=searchHits, verbosity=0, 
                      rearrangements="NNI", suboptimal=1e-08, ...) {
  if (is.null(treeOrder <- attr(tree, 'order')) || treeOrder != 'preorder') tree <- Preorder(tree)
   
  epsilon <- 1e-08
  if (is.null(attr(tree, "score"))) attr(tree, "score") <- TreeScorer(tree, dataset, ...)
  bestScore <- attr(tree, "score")
  if (verbosity >= 0) cat("\n* Initial score:", bestScore)
  if (returnAll) {
    nullForest <- vector('list', ratchIter)
    forest <- nullForest
    forestScores <- rep(NA, ratchIter)
  }

  bestScoreHits <- 0
  iterationsCompleted <- 0
  for (i in 1:ratchIter) {
    if (verbosity >= 1) cat ("\n - Running NNI on bootstrapped dataset. ")
    bootstrapTree <- BootstrapTree(tree, dataset, TreeScorer=TreeScorer,
                                   maxIter=bootstrapIter, maxHits=bootstrapHits, 
                                   rooted=rooted, verbosity=verbosity-2, ...)
    
    if (verbosity >= 1) cat ("\n - Rearranging from new candidate tree:")
    if (is.character(rearrangements)) {
      Rearrangements <- if (rooted) {
        if (rearrangements == "TBR") 
          list(RootedTBR, RootedSPR, RootedNNI) else 
        if (rearrangements == "TBR only") 
          list(RootedTBR) else
        if (rearrangements == "SPR") 
          list(RootedSPR, RootedNNI) else 
        if (rearrangements == "SPR only") 
          list(RootedSPR) else list(RootedNNI)
      } else {
        if (rearrangements == "TBR") 
          list(TBR, SPR, NNI) else 
        if (rearrangements == "TBR only") 
          list(TBR) else
        if (rearrangements == "SPR") 
          list(SPR, NNI) else 
        if (rearrangements == "SPR only") 
          list(SPR) else list(NNI)
      }
    } else Rearrangements <- rearrangements
    
    candidate <- bootstrapTree
    for (Func in Rearrangements) {
      candidate <- TreeSearch(candidate, dataset, TreeScorer=TreeScorer, Rearrange=Func, 
                                verbosity=verbosity, maxIter=searchIter, maxHits=searchHits, ...)
    }
    
    candScore <- attr(candidate, 'score')
    if (verbosity >= 1) cat("\n - Candidate tree scored ", candScore)
    if (returnAll && candScore < (bestScore + suboptimal)) { # Worth saving this tree in forest
      forest[[i]] <- candidate
      forestScores[i] <- candScore
    }
    if ((candScore + epsilon) < bestScore) {
      # New 'best' tree
      tree <- candidate
      bestScore <- candScore
      bestScoreHits <- 1L
      
    } else if (bestScore + epsilon > candScore) { # i.e. best == cand, allowing for floating point error
      bestScoreHits <- bestScoreHits + 1
      tree <- candidate
    }
    if (verbosity >= 0) cat("\n* Best score after", i, "/", ratchIter, "ratchet iterations:", bestScore, "( hit", bestScoreHits, "/", ratchHits, ")")
    if (bestScoreHits >= ratchHits) {
      iterationsCompleted <- i
      break()
    }
  } # end for
  if (iterationsCompleted == 0) iterationsCompleted <- ratchIter
  if (verbosity >= 0) cat ("\nCompleted parsimony ratchet after", iterationsCompleted, "iterations with score", bestScore, "\n")
   
  if (returnAll) {
    keepers <- !is.na(forestScores) & forestScores < bestScore + suboptimal
    forestScores <- forestScores[keepers]
    forest <- forest[keepers]
    if (verbosity >=0 ) cat("\n - Keeping", sum(keepers), "trees from iterations numbered:\n   ", which(keepers))
    if (length(forest) > 1) {
      class(forest) <- 'multiPhylo'
      ret <- unique(forest)
      if (verbosity >=0) cat("\n - Removing duplicates leaves", length(ret), "unique trees")
    } else if (length(forest) == 1) {
      class(forest) <- 'phylo'
      ret <- forest
    } else {
      stop('\nNo trees!? Is suboptimal set to a sensible (positive) value?')
    }
    scores.unique <- vapply(ret, attr, double(1), 'score')
    cat('\nFound', sum(scores.unique == min(scores.unique)), 'unique MPTs and', length(ret) - sum(scores.unique == min(scores.unique)), 'suboptimal trees.\n')
  } else {
    ret <- tree
    attr(ret, 'hits') <- NULL
  }
  return (ret)
}