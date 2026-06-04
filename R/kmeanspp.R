#' k-means++ clustering
#'
#' k-means++ clustering \insertCite{Arthur2007}{TreeDist} improves the speed and
#' accuracy of standard \code{\link[stats]{kmeans}} clustering
#' \insertCite{Hartigan1979}{TreeDist} by preferring initial cluster centres 
#' that are far from others.
#' A scalable version of the algorithm has been proposed for larger data sets
#' \insertCite{Bahmani2012}{TreeDist}, but is not implemented here.
#' 
#' @param x Numeric matrix of data, or an object that can be coerced to such a
#' matrix (such as a numeric vector or a data frame with all numeric columns).
#' @param k Integer specifying the number of clusters, _k_.
#' @param nstart Positive integer specifying how many random sets should be
#' chosen
#' @param \dots additional arguments passed to \code{\link[stats]{kmeans}}
#' @references
#' \insertAllCited{} 
#' 
#' @template MRS
#' @seealso \code{\link[stats]{kmeans}}
#' @examples
#' # Generate random points
#' set.seed(1)
#' x <- cbind(c(rnorm(10, -5), rnorm(5, 1), rnorm(10, 6)),
#'            c(rnorm(5, 0), rnorm(15, 4), rnorm(5, 0)))
#' 
#' # Conventional k-means may perform poorly
#' klusters <- kmeans(x, cent = 5)
#' plot(x, col = klusters$cluster, pch = rep(15:19, each = 5))
#' 
#' # Here, k-means++ recovers a better clustering
#' plusters <- KMeansPP(x, k = 5)
#' plot(x, col = plusters$cluster, pch = rep(15:19, each = 5))
#' @family cluster functions
#' @importFrom stats kmeans
#' @export
KMeansPP <- function(x, k = 2, nstart = 10, ...) UseMethod("KMeansPP")

#' @export
KMeansPP.matrix <- function(x, k = 2, nstart = 10, ...) {
  if (k < 2) {
    return(kmeans(x, centers = k, ...))
  }
  
  n <- dim(x)[[1]]
  nCol <- dim(x)[[2]]
  ret <- list(tot.withinss = Inf)

  # k-means++ seeding needs, at each step, only the Euclidean distance from the
  # chosen centre to all n points: an O(n * nCol) computation per centre.
  # Computing these rows on the fly avoids materializing the full n * n distance
  # matrix (O(n ^ 2) time and memory), so this method scales to large n.
  .DistanceRow <- function(centre) {
    sqrt(.rowSums((x - rep(x[centre, ], each = n)) ^ 2, n, nCol))
  }

  for (start in seq_len(nstart)) {
    centres <- integer(k)
    centres[1L] <- sample.int(n, 1L)
    min_d <- .DistanceRow(centres[1L])

    for (i in 2:k) {
      p <- min_d ^ 2
      if (!any(p != 0)) {
        stop("Not enough distinct data points to compute clustering")
      }
      centres[i] <- sample.int(n, 1L, prob = p)
      min_d <- pmin.int(min_d, .DistanceRow(centres[i]))
    }

    proposal <- kmeans(x, centers = x[centres, ], ...)
    if (proposal[["tot.withinss"]] < ret[["tot.withinss"]]){
      ret <- proposal
    }
  }
  
  # Return:
  ret
}

#' @export
KMeansPP.numeric <- function(x, k = 2, nstart = 10, ...) {
  KMeansPP(matrix(x, ncol = 1), k = k, nstart = nstart, ...)
}

#' @export
KMeansPP.dist <- function(x, k = 2, nstart = 10, ...) {
  if (k < 2) {
    return(kmeans(x, centers = k, ...))
  }
  
  n <- attr(x, "Size")
  ret <- list(tot.withinss = Inf)
  # This method clusters the n rows of the n * n distance matrix as
  # n-dimensional vectors, so it is O(n ^ 2) by construction.  Callers that hold
  # point coordinates should use KMeansPP.matrix(), which seeds in O(n * dim).
  d <- as.matrix(x)

  for (start in seq_len(nstart)) {
    centres <- integer(k)
    centres[1L] <- sample.int(n, 1L)
    min_d <- d[centres[1L], ]

    for (i in 2:k) {
      p <- min_d ^ 2
      if (!any(p != 0)) {
        stop("Not enough distinct data points to compute clustering")
      }
      centres[i] <- sample.int(n, 1L, prob = p)
      min_d <- pmin.int(min_d, d[centres[i], ])
    }

    # Pass the already-materialized matrix d: kmeans() coerces its data argument
    # with as.matrix(), so passing the dist object x would build a second
    # n * n matrix needlessly.
    proposal <- kmeans(d, centers = d[centres, ], ...)
    if (proposal[["tot.withinss"]] < ret[["tot.withinss"]]){
      ret <- proposal
    }
  }
  
  # Return:
  ret
}
