
#' @keywords internal
#' @export
Min <- function (x, inappLevel) {
  if (length(inappLevel)) return(sum(2^(c(0:(inappLevel-2), inappLevel:12)) %in% unique(x)))
  return (sum(2^(0:12) %in% unique(x)))
}


#' Quick sample
#' 
#' Faster than inbuilt sample because it avoids some checks
#' @param x vector to sample
#' @param len length of vector
#' @keywords internal
#' @export
SampleOne <- function (x, len = length(x)) x[sample.int(len, 1L, FALSE, NULL, FALSE)]
