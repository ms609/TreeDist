#' k-means++ clustering
#'
#' k-means++ clustering \insertCite{Arthur2007}{TreeDist} improves the speed and
#' accuracy of standard \code{\link[stats]{kmeans}} clustering by selecting 
#' random seeds that are far from existing cluster centroids.
#' 
#' @param x Numeric matrix of data, or an object that can be coerced to such a
#' matrix (such as a numeric vector or a data frame with all numeric columns).
#' @param k Integer specifying the number of clusters, _k_.
#' @param nstart Positive integer specifying how many random sets should be
#' chosen
#' @param \dots additional arguments passed to \code{\link[stats]{kmeans}}
#' @references
#' \insertAllCited{}#'  
#'  
#' @author Modified by Martin R. Smith from [code](
#' https://rdrr.io/cran/LICORS/src/R/kmeanspp.R) by Georg M. Goerg.
#' @export
#' @seealso \code{\link[stats]{kmeans}}
#' @examples 
#' # Generate random points
#' set.seed(1)
#' n <- 100
#' x <- rbind(
#'   matrix(rnorm(n), ncol = 2),
#'   matrix(runif(n * 2, -1, 1), ncol = 2)
#' )
#' 
#' # Compute clustering
#' clusters <- kmeanspp(x, k = 5)
#' 
#' # Plot
#' plot(x, col = clusters$cluster, pch = 16)
#' @importFrom stats kmeans
kmeanspp <- function(x, k = 2, nstart = 10, ...) {
  
  if (length(dim(x)) == 0) {
    x <- matrix(x, ncol = 1)
  } else {
    x <- cbind(x)
  }
  
  if (k < 2) {
    return(kmeans(x, centers = k, ...))
  }
  
  n <- dim(x)[1]
  ret <- list(tot.withinss = Inf)
  
  for (restart in seq_len(nstart)) {  
    centres <- numeric(k)
    centres[1:2] <- sample.int(n, 1)
    
    for (i in seq_len(k)[-1]) {
      dists <- apply(x[centres, , drop = FALSE], 1, 
                     function(centre) rowSums((x - centre) ^ 2))
      probs <- apply(dists, 1, min)
      probs[centres] <- 0
      centres[i] <- sample.int(n, 1, prob = probs)
    }
    
    proposal <- kmeans(x, centers = x[centres, ], ...)
    if (proposal$tot.withinss < ret$tot.withinss){
      ret <- proposal
    }
  }
  
  # Return:
  ret
}
