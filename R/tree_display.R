#' Consensus without taxa
#' 
#' Displays a consensus plot with selected taxa excluded.
#' 
#' A useful way to gain resolution if a few wildcard taxa obscure a consistent
#' set of relationship.
#' 
#' @param trees A list of phylogenetic trees, of class `multiPhylo` or `list`
#' @param tip A character vector specifying the names (or numbers) of tips to
#'                drop (using ape::drop.tip)
#' @param \dots Additional parameters to pass to ape::[consensus]
#'                
#' @return A consensus tree without the excluded taxa
#' @author Martin R. Smith
#' @importFrom ape consensus drop.tip
#' @export
ConsensusWithout <- function (trees, tip, ...) {
  if (class(trees) == 'phylo') {
    drop.tip(trees, tip=tip) 
  } else {
    consensus(lapply(trees, drop.tip, tip=tip), ...)
  }
}

#' @describeIn ConsensusWithout Adds missing taxa to a plotted consensus tree
#' @param position Where to plot the missing taxa.  See [legend] for options.
#' @param \dots Other parameters to pass to [legend]
#' @importFrom graphics legend
#' @export
#' @author Martin R. Smith
MarkMissing <- function (tip, position='bottomleft', ...) {
  legend(position, legend=gsub('_', ' ', tip, fixed=TRUE),
         lwd=1, lty=2, bty='n', ...)
}