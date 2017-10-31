#' @param mergeEdges (optional) vector of length 1 or 2, listing edge(s) to be joined:
#'                   In SPR, this is where the pruned subtree will be reconnected.
#'                   In TBR, these edges will be reconnected (so must be on opposite
#'                   sides of \code{edgeToBreak}); if only a single edge is specified,
#'                   the second will be chosen at random
