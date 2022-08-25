#' Cluster size statistics
#' @name cluster-statistics
#' 
#' @param x Matrix in which each row lists the coordinates of a point
#' in a Euclidian space; or, where supported, `dist` object specifying
#' distances between each pair of points.
#' @param cluster Optional integer vector specifying the cluster or group to
#' which each row in `x` belongs.
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
SumOfRanges <- function(x, cluster = 1) {
  if (is.null(dim(x))) {
    warning(paste0("`x` lacks dimensions. ",
                   "Did you subset without specifying `drop = FALSE`?"))
    x <- matrix(x, 1)
  }
  
  # Return:
  vapply(seq_along(unique(cluster)), 
         function(i) .SumOfRanges(x[cluster == i, , drop = FALSE]),
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
#' @importFrom stats var
#' @export
SumOfVariances <- function(x, cluster = 1) {
  if (is.null(dim(x))) {
    warning(paste0("`x` lacks dimensions. ",
                   "Did you subset without specifying `drop = FALSE`?"))
    x <- matrix(x, 1)
  }
  
  # Return:
  vapply(seq_along(unique(cluster)), 
         function(i) .SumOfVariances(x[cluster == i, , drop = FALSE]),
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
MeanCentroidDistance <- function(x, cluster = 1) {
  if (is.null(dim(x))) {
    warning(paste0("`x` lacks dimensions. ",
                   "Did you subset without specifying `drop = FALSE`?"))
    x <- matrix(x, 1)
  }
  
  # Return:
  vapply(seq_along(unique(cluster)),
         function(i) .MeanCentroidDist(x[cluster == i, , drop = FALSE]),
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
#' @return `DistanceFromMedian()` returns a numeric specifying the mean distance
#' of each point (except the median) from the median point of its cluster.
#' @examples DistanceFromMedian(points, cluster)
#' @export
DistanceFromMedian <- function(x, cluster = 1) UseMethod("DistanceFromMedian")

#' @rdname cluster-statistics
#' @export
DistFromMed <- DistanceFromMedian

#' @export
DistanceFromMedian.dist <- function(x, cluster = 1) {
  d <- as.matrix(x)
  vapply(seq_along(unique(cluster)),
         function(i) {
           .DistanceFromMedian.dist(
             d[cluster == i, cluster == i, drop = FALSE])
         },
         numeric(1))
}

#' @export
DistanceFromMedian.numeric <- function(x, cluster = 1) {
  if (is.null(dim(x))) {
    warning(paste0("`x` lacks dimensions. ",
                   "Did you subset without specifying `drop = FALSE`?"))
    x <- matrix(x, 1)
  }
  
  # Return:
  vapply(seq_along(unique(cluster)),
         function(i) .DistanceFromMedian(x[cluster == i, , drop = FALSE]),
         numeric(1))
}

.DistanceFromMedian <- function(x) {
  if (dim(x)[1] > 1) {
    .DistanceFromMedian.dist(as.matrix(dist(x)))
  } else {
    NA_real_
  }
}

.DistanceFromMedian.dist <- function(d) {
  if (dim(d)[1] > 1) {
    medPoint <- which.min(unname(colSums(d)))
    mean(d[medPoint, -medPoint])
  } else {
    NA_real_
  }
}

#' @rdname cluster-statistics
#' @return `MeanNN()` returns a numeric specifying the mean distance from each
#' point within a cluster to its nearest neighbour.
#' @examples MeanNN(points, cluster)
#' @export
MeanNN <- function(x, cluster = 1) UseMethod("MeanNN")

#' @export
MeanNN.dist <- function(x, cluster = 1) {
  d <- as.matrix(x)
  diag(d) <- NA_real_
  vapply(seq_along(unique(cluster)),
         function(i) .MeanNN.dist(d[cluster == i, cluster == i, drop = FALSE]),
         numeric(1))
}

.MeanNN.dist <- function(x) {
  if (dim(x)[1] > 1) {
    mean(apply(x, 1, min, na.rm = TRUE))
  } else {
    NA_real_
  }
}

#' @export
MeanNN.numeric <- function(x, cluster = 1) {
  if (is.null(dim(x))) {
    warning(paste0("`x` lacks dimensions. ",
                   "Did you subset without specifying `drop = FALSE`?"))
    x <- matrix(x, 1)
  }
  
  # Return:
  vapply(seq_along(unique(cluster)),
         function(i) .MeanNN(x[cluster == i, , drop = FALSE]),
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
MeanMSTEdge <- function(x, cluster = 1) UseMethod("MeanMSTEdge")

#' @export
MeanMSTEdge.dist <- function(x, cluster = 1) {
  d <- as.matrix(x)
  diag(d) <- NA_real_
  vapply(seq_along(unique(cluster)),
         function(i) {
           .MeanMSTEdge.dist(d[cluster == i, cluster == i, drop = FALSE])
         },
         numeric(1))
}

.MeanMSTEdge.dist <- function(x) {
  n <- dim(x)[1]
  # Return:
  if (n > 1) {
    MSTLength(as.dist(x)) / (n - 1)
  } else {
    NA_real_
  }
}

#' @export
MeanMSTEdge.numeric <- function(x, cluster = 1) {
  if (is.null(dim(x))) {
    warning(paste0("`x` lacks dimensions. ",
                   "Did you subset without specifying `drop = FALSE`?"))
    x <- matrix(x, 1)
  }
  
  # Return:
  vapply(seq_along(unique(cluster)),
         function(i) .MeanMSTEdge(x[cluster == i, , drop = FALSE]),
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
