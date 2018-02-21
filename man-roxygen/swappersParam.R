#' @param swappers A list of functions to use to conduct edge rearrangement during tree search.
#'                 Provide functions like \code{\link{NNISwap}} to shuffle root position,
#'                 or \code{\link{RootedTBRSwap}} if the position of the root should be retained.
#'                 You may wish to use extreme swappers (such as \acronym{TBR}) early in the list,
#'                 and a more subtle rearranger (such as \acronym{NNI}) later in the list to make
#'                 incremental tinkerings once an almost-optimal tree has been found.
