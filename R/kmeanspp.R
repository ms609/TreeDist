#' k-means++ clustering
#'
#' k-means++ clustering \insertCite{Arthur2007}{TreeDist} improves the speed and
#' accuracy of standard \code{\link[stats]{kmeans}} clustering by preferring
#' initial cluster centres that are far from others.
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
  
  n <- dim(x)[1]
  ret <- list(tot.withinss = Inf)
  d <- as.matrix(dist(x))
  
  for (start in seq_len(nstart)) {  
    centres <- c(sample.int(n, 1), numeric(k - 1))
    
    for (i in 2:k) {
      centres[i] <- sample.int(
        n, 1, prob = apply(d[centres, , drop = FALSE], 2, min) ^ 2
      )
    }
    
    proposal <- kmeans(x, centers = x[centres, ], ...)
    if (proposal$tot.withinss < ret$tot.withinss){
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
  d <- as.matrix(x)
  
  for (start in seq_len(nstart)) {  
    centres <- c(sample.int(n, 1), numeric(k - 1))
    
    for (i in 2:k) {
      centres[i] <- sample.int(
        n, 1, prob = apply(d[centres, , drop = FALSE], 2, min) ^ 2
      )
    }
    
    proposal <- kmeans(x, centers = d[centres, ], ...)
    if (proposal$tot.withinss < ret$tot.withinss){
      ret <- proposal
    }
  }
  
  # Return:
  ret
}
