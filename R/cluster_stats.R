#' Cluster size statistics
#' @name cluster-statistics
#' 
#' @param mapping Matrix in which each row lists the coordinates of a point
#' in a Euclidian space.
#' @param cluster Optional integer vector specifying the cluster or group to
#' which each row in `mapping` belongs.
#' 
#' @family tree space functions
#' @template MRS
NULL

#' @return `SumOfRanges()` returns a numeric specifying the sum of ranges
#' within each cluster across all dimensions.
#' @rdname cluster-statistics
#' @export
SumOfRanges <- function(mapping, cluster = 1) {
  if (is.null(dim(mapping))) {
    warning(paste0("`mapping` lacks dimensions. ",
                   "Did you subset without specifying `drop = FALSE`?"))
    mapping <- matrix(mapping, 1)
  }
  
  # Return:
  vapply(seq_along(unique(cluster)), 
         function(i) .SumOfRanges(mapping[cluster == i, , drop = FALSE]),
         numeric(1))
}

.SumOfRanges <- function(x) {
  ranges <- apply(x, 2, range)
  
  # Return:
  sum(ranges[2, ] - ranges[1, ])
}

#' @rdname cluster-statistics
#' @return `SumOfRanges()` returns a numeric specifying the sum of variances
#' within each cluster across all dimensions.
#' @export
SumOfVariances <- function(mapping, cluster = 1) {
  if (is.null(dim(mapping))) {
    warning(paste0("`mapping` lacks dimensions. ",
                   "Did you subset without specifying `drop = FALSE`?"))
    mapping <- matrix(mapping, 1)
  }
  
  # Return:
  vapply(seq_along(unique(cluster)), 
         function(i) .SumOfVariances(mapping[cluster == i, , drop = FALSE]),
         numeric(1))
}

#' @rdname cluster-statistics
#' @export
SumOfVars <- SumOfVariances

.SumOfVariances <- function(x) {
  sum(apply(x, 2, var))
}


.Apply <- if (packageVersion("base") < "4.1.0") {
  function(x, ...) {
    if (dim(x)) {
      apply(x, ...)
    } else {
      matrix(apply(x, ...), 1)
    }
  }
} else {
  function(x, ...) apply(x, ..., simplify = FALSE)
}
