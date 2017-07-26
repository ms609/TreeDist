#' Parsimony Ratchet
#'
#' \code{Ratchet} uses the parsimony ratchet (Nixon 1999) to search for a more parsimonious tree.
#'
#' @template treeParam 
#' @param data a dataset in the format required by TreeScorer
#' @template concavityParam
#' @param all Set to \code{TRUE} to report all MPTs encountered during the search, perhaps to analyze consensus
#' @param outgroup a vector specifying all tips in the outgroup; if unspecified then identical trees with different roots will be considered unique;
#' @param maxit   maximum ratchet iterations to perform;
#' @param maxIter maximum rearrangements to perform on each bootstrap or ratchet iteration;
#' @param maxHits maximum times to hit best score before terminating a tree search within a pratchet iteration;
#' @param k stop when k ratchet iterations have found the same best score;
#' @param verbosity larger numbers provides more verbose feedback to the user;
#' @param rearrangements method(s) to use when rearranging trees: 
#'        a vector containing any combination of the strings "NNI", "SPR" or "TBR";
#' @param \dots other arguments to pass to subsequent functions.
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
#' @seealso \code{\link{SectorialSearch}}
#' 
#' @examples Ratchet(RandomTree(Lobo.phy), SigSut.phy, outgroup='Cricocosmia')
#' 
#' @keywords  tree 
#' @export
Ratchet <- function (tree, data, TreeScorer=FitchScore, all=FALSE, outgroup=NULL, 
                      pratchiter=100, searchiter=5000, searchhits=40, pratchhits=10, track=0, 
                      rearrangements="NNI", suboptimal=1e-08, ...) {
  epsilon <- 1e-08
  if (is.null(attr(tree, "score"))) attr(tree, "score") <- TreeScorer(tree, data)
  best.score <- attr(tree, "score")
  if (track >= 0) cat("\n* Initial score:", best.score)
  if (all) {
    null.forest <- vector('list', pratchiter)
    forest <- null.forest
    forest.scores <- rep(NA, pratchiter)
  }

  best.score.hits <- 0
  iterations.completed <- 0
  for (i in 1:pratchiter) {
    if (track >= 0) cat ("\n - Running NNI on bootstrapped dataset. ")
    bstree <- BootstrapTree(phy=tree, x=data, maxIter=searchiter, maxHits=searchhits,
                        TreeScorer=TreeScorer,  track=track - 1, ...)
    
    if (track >= 0) cat ("\n - Running", ifelse(is.null(rearrangements), "NNI", rearrangements), "from new candidate tree:")
    if (rearrangements == "TBR") {
      candidate <- DoTreeSearch(bstree,    data, TreeScorer=TreeScorer, method='TBR', track=track, maxIter=searchiter, maxHits=searchhits, ...)
      candidate <- DoTreeSearch(candidate, data, TreeScorer=TreeScorer, method='SPR', track=track, maxIter=searchiter, maxHits=searchhits, ...)
      candidate <- DoTreeSearch(candidate, data, TreeScorer=TreeScorer, method='NNI', track=track, maxIter=searchiter, maxHits=searchhits, ...)
    } else if (rearrangements == "TBR only") {  
      candidate <- DoTreeSearch(bstree,    data, TreeScorer=TreeScorer, method='TBR', track=track, maxIter=searchiter, maxHits=searchhits, ...)
    } else if (rearrangements == "SPR") {       
      candidate <- DoTreeSearch(bstree,    data, TreeScorer=TreeScorer, method='SPR', track=track, maxIter=searchiter, maxHits=searchhits, ...)
      candidate <- DoTreeSearch(candidate, data, TreeScorer=TreeScorer, method='NNI', track=track, maxIter=searchiter, maxHits=searchhits, ...)
    } else if (rearrangements == "SPR only") {  
      candidate <- DoTreeSearch(bstree,    data, TreeScorer=TreeScorer, method='SPR', track=track, maxIter=searchiter, maxHits=searchhits, ...)
    } else {  
      candidate <- DoTreeSearch(bstree,    data, TreeScorer=TreeScorer, method='NNI', track=track, maxIter=searchiter, maxHits=searchhits, ...)
    }
    cand.score <- attr(candidate, 'score')
    if ((cand.score + epsilon) < best.score) {
      # New 'best' tree
      if (all) {
        forest[[i]] <- if (is.null(outgroup)) candidate else Root(candidate, outgroup)
        forest.scores[i] <- cand.score
      }
      tree <- candidate
      best.score <- cand.score
      best.score.hits <- 1
    } else if (best.score + epsilon > cand.score) { # i.e. best == cand, allowing for floating point error
      best.score.hits <- best.score.hits + 1
      tree <- candidate
      if (all) {
        forest[[i]] <- if (is.null(outgroup)) candidate else Root(candidate, outgroup)
        forest.scores[i] <- cand.score
      }
    } else if (cand.score < (best.score + suboptimal) && all) {
      forest[[i]] <- if (is.null(outgroup)) candidate else Root(candidate, outgroup)
      forest.scores[i] <- cand.score
    }
    if (track >= 0) cat("\n* Best score after", i, "/", pratchiter, "pratchet iterations:", best.score, "( hit", best.score.hits, "/", pratchhits, ")")
    if (best.score.hits >= pratchhits) {
      iterations.completed <- i
      break()
    }
  } # end for
  if (iterations.completed == 0) iterations.completed <- pratchiter
  if (track >= 0) cat ("\nCompleted parsimony ratchet after", iterations.completed, "iterations with score", best.score, "\n")
   
  if (all) {
    keepers <- !is.na(forest.scores) & forest.scores < best.score + suboptimal
    forest.scores <- forest.scores[keepers]
    forest <- forest[keepers]
    if (length(forest) > 1) {
      class(forest) <- 'multiPhylo'
      ret <- unique(forest)
    } else if (length(forest) == 1) {
      class(forest) <- 'phylo'
      ret <- forest
    } else {
      stop('No trees!? Is suboptimal set to a sensible (positive) value?')
    }
    scores.unique <- vapply(ret, attr, double(1), 'score')
    cat('Found', sum(scores.unique == min(scores.unique)), 'unique MPTs and', length(ret) - sum(scores.unique == min(scores.unique)), 'suboptimal trees.\n')
    if (is.null(outgroup)) warning('"outgroup" not specified, so some "unique" trees may have same topology but distinct roots.')
  } else {
    ret <- tree
    attr(ret, 'hits') <- NULL
  }
  return (ret)
}