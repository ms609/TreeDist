
#' @keywords internal
#' @export
Min <- function (x, inappLevel) {
  if (length(inappLevel)) return(sum(2^(c(0:(inappLevel-2), inappLevel:12)) %in% unique(x)))
  return (sum(2^(0:12) %in% unique(x)))
}