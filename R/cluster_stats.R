#' Cluster size statistics
#' @name cluster-statistics
#' 
#' @param mapping Matrix in which each row lists the coordinates of a point
#' in a Euclidian space.
#' @param cluster Optional integer vector specifying the cluster or group to
#' which each row in `mapping` belongs.
#' 
#' @examples
#' points <- rbind(matrix(1:16, 4), rep(1, 4), matrix(1:32, 8, 4) / 10)
#' cluster <- rep(1:3, c(4, 1, 8))
#' 
#' plot(
#'   points[, 1:2], # Plot first two dimensions of four-dimensional space
#'   col = cluster, pch = cluster, # Style by cluster membership
#'   asp = 1, # Fix aspect ratio to avoid distortion
#'   ann = FALSE, frame = FALSE # Simple axes
#' )
#' 
#' @family tree space functions
#' @template MRS
NULL

#' @return `SumOfRanges()` returns a numeric specifying the sum of ranges
#' within each cluster across all dimensions.
#' @rdname cluster-statistics
#' @examples SumOfRanges(points, cluster)
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
#' @return `SumOfVariances()` returns a numeric specifying the sum of variances
#' within each cluster across all dimensions.
#' @examples SumOfVariances(points, cluster)
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


#' @rdname cluster-statistics
#' @return `MeanCentroidDistance()` returns a numeric specifying the mean
#' distance from the centroid to points in each cluster.
#' @examples MeanCentroidDistance(points, cluster)
#' @export
MeanCentroidDistance <- function(mapping, cluster = 1) {
  if (is.null(dim(mapping))) {
    warning(paste0("`mapping` lacks dimensions. ",
                   "Did you subset without specifying `drop = FALSE`?"))
    mapping <- matrix(mapping, 1)
  }
  
  # Return:
  vapply(seq_along(unique(cluster)),
         function(i) .MeanCentroidDist(mapping[cluster == i, , drop = FALSE]),
         numeric(1))
}

#' @rdname cluster-statistics
#' @export
MeanCentDist <- MeanCentroidDistance

#' @rdname cluster-statistics
#' @export
MeanCentroidDist <- MeanCentroidDistance

.MeanCentroidDist <- function(x) {
  recentred <- t(t(x) - apply(x, 2, mean))
  
  # Return:
  mean(sqrt(rowSums(recentred ^ 2)))
}

#' @rdname cluster-statistics
#' @return `MeanNN()` returns a numeric specifying the mean distance from each
#' point within a cluster to its nearest neighbour.
#' @examples MeanNN(points, cluster)
#' @export
MeanNN <- function(mapping, cluster = 1) {
  if (is.null(dim(mapping))) {
    warning(paste0("`mapping` lacks dimensions. ",
                   "Did you subset without specifying `drop = FALSE`?"))
    mapping <- matrix(mapping, 1)
  }
  
  # Return:
  vapply(seq_along(unique(cluster)),
         function(i) .MeanNN(mapping[cluster == i, , drop = FALSE]),
         numeric(1))
}

.MeanNN <- function(x) {
  if (dim(x)[1] > 1) {
    distances <- as.matrix(dist(x))
    diag(distances) <- NA_real_
    
    # Return:
    mean(apply(distances, 1, min, na.rm = TRUE))
    
  } else {
    # Return:
    NA_real_
  }
}

#' @rdname cluster-statistics
#' @return `MeanMSTEdge()` returns a numeric specifying the mean length of an
#' edge in the minimum spanning tree of points within each cluster.
#' @examples MeanMSTEdge(points, cluster)
#' @export
MeanMSTEdge <- function(mapping, cluster = 1) {
  if (is.null(dim(mapping))) {
    warning(paste0("`mapping` lacks dimensions. ",
                   "Did you subset without specifying `drop = FALSE`?"))
    mapping <- matrix(mapping, 1)
  }
  
  # Return:
  vapply(seq_along(unique(cluster)),
         function(i) .MeanMSTEdge(mapping[cluster == i, , drop = FALSE]),
         numeric(1))
}

#' @importFrom TreeTools MSTLength
.MeanMSTEdge <- function(x) {
  n <- dim(x)[1]
  # Return:
  if (n > 1) {
    MSTLength(dist(x)) / (n - 1)
  } else {
    NA_real_
  }
}
