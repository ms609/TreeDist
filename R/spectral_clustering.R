#' Eigenvalues for spectral clustering
#'
#' Spectral clustering emphasizes nearest neighbours when forming clusters;
#' it avoids some of the issues that arise from clustering around means / 
#' medoids.
#'
#' @param D Square matrix or `dist` object containing Euclidean distances
#'   between data points.
#' @param nn Integer specifying number of nearest neighbours to consider
#' @param nEig Integer specifying number of eigenvectors to retain.
#' @author Adapted by MRS from script by [Nura
#'   Kawa](https://rpubs.com/nurakawa/spectral-clustering)
#' @return `SpectralEigens()` returns spectral eigenvalues that can then be
#'   clustered using a method of choice.
#' @examples
#' library('TreeTools', quietly = TRUE, warn.conflict = FALSE)
#' trees <- as.phylo(0:18, nTip = 8)
#' distances <- ClusteringInfoDistance(trees)
#' eigens <- SpectralEigens(distances)
#' # Perform clustering:
#' clusts <- kmeans(dist(eigens), centers = 3)
#' plot(eigens, pch = 15, col = clusts$cluster)
#' plot(cmdscale(distances), pch = 15, col = clusts$cluster)
#' @family tree space functions
#' @export
SpectralEigens <- function (D, nn = 10L, nEig = 2L) {
  
  MutualKnnGraph <- function (D, nn) {
    D <- as.matrix(D)
    dims <- dim(D)
    
    # intialize the knn matrix
    knn_mat <- matrix(FALSE, nrow = dims[1], ncol = dims[2])
    
    # find the 10 nearest neighbors for each point
    for (i in seq_len(nrow(D))) {
      neighbor_index <- order(D[i, ])[2:(nn + 1)]
      knn_mat[i, ][neighbor_index] <- TRUE
    }
    
    # Now we note that i,j are neighbors iff K[i,j] = 1 or K[j,i] = 1
    knn_mat <- knn_mat | t(knn_mat) # find mutual knn
    
    # Return:
    knn_mat
  }
  
  GraphLaplacian <- function(W) {
    stopifnot(nrow(W) == ncol(W))
    
    g <- colSums(W) # degrees of vertices
    n <- dim(W)[1]
    D_half <- diag(1 / sqrt(g))
    
    # Return:
    diag(n) - D_half %*% W %*% D_half
  }
  
  W <- MutualKnnGraph(D, nn) # 1. matrix of similarities
  L <- GraphLaplacian(W) # 2. compute graph laplacian
  ei <- eigen(L, symmetric = TRUE) # 3. Compute the eigenvectors and values of L
  
  # Return the eigenvectors of the n_eig smallest eigenvalues:
  ei$vectors[, nrow(L) - rev(seq_len(nEig))]
}

#' @export
#' @rdname SpectralEigens
SpectralClustering <- function (D, nn = 10L, nEig = 2L) {
  .Deprecated("SpectralEigens")
  SpectralEigens(D, nn, nEig)
}