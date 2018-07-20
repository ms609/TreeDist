#' @param stopAtPeak Logical specifying whether to terminate search once a subsequent 
#'                   iteration recovers a sub-optimal score. Useful with methods that
#'                   return all trees one rearrangement from the current tree, such 
#'                   as \code{\link{AllTBR}}.  Will be overridden if a passed function
#'                   has an attribute `stopAtPeak` set by
#'                   `attr(FunctionName, 'stopAtPeak') <- TRUE`.
