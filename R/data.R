#' One Overlap
#' 
#' Information content of a single overlap.
#' 
#' Internal cache of all values for up to 100 leaves, called from within 
#' `SharedPhylogeneticInfoSplits`.
#' 
#' The benefits of caching decrease as `nLeaves` increases, 
#' as a match is less likely; moreover, the cost, in terms of storage
#' space, rises.
#' 
#' @noRd
#' @keywords datasets
"oneOverlap"