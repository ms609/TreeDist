#' Bootstrap tree search with inapplicable data
#' 
#' @template labelledTreeParam

#'
#' @return A tree that is optimal under a random sampling of the original characters
#' @export  
BootstrapTree <- function (phy, x, maxiter, maxhits, ParsimonyScorer = phangorn::fitch, track=1, ...) {
## Simplified version of phangorn::bootstrap.phyDat, with bs=1 and multicore=FALSE
  at <- attributes(x)
  bootstrapped <- BootstrapWeightings(at)
  keep <- bootstrapped > 0
  x <- x[keep, ] ## TODO use TipsAreColumns family
  attr(x, 'weight') <- bootstrapped[which(keep)]
  attr(x, 'min.steps') <- at$min.steps[keep]
  attr(x, 'info.amounts') <- at$info.amounts[keep, ]
  attr(x, 'unique.tokens') <- at$unique.tokens[keep]
  attr(x, 'sa.weights') <- at$sa.weights[keep]
  attr(x, 'nr') <- sum(keep)
  attr(x, 'inapp.level') <- at$inapp.level
  attr(phy, 'score') <- NULL
  class(x) <- 'fitchDat'
  res <- DoTreeSearch(phy, x, ParsimonyScorer=ParsimonyScorer, method='NNI', maxiter=maxiter,
                    maxhits=maxhits, track=track-1, ...)
  attr(res, 'score') <- NULL
  attr(res, 'hits') <- NULL
  res
}

#' @keywords internal
#' @export
BootstrapWeightings <- function (attribs) {
  weight <- attribs$weight
  v <- rep(1:length(weight), weight)
  BS <- tabulate(sample(v, replace=TRUE), length(weight)) 
}