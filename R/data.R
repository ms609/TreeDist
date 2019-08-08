#' One Overlap
#' 
#' Information content of a single overlap
#' Internal cache of all values for up to 100 terminals, called from within `MutualPhylogeneticInfoSplits`.
#' 
#' The benefits of cacheing decrease an nTerminals increases (as a match is less likely)
#' and the cost rises (in terms of storage space).
#' 
#' @keywords datasets
"oneOverlap"