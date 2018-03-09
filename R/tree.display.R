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
#'                
#'  @return A consensus tree without the excluded taxa
#'  @author Martin R. Smith
#'  @importFrom ape consensus drop.tip
#'  @export
ConsensusWithout <- function (trees, tip) {
  consensus(lapply(trees, drop.tip, tip=tip))
}

