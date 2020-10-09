#' Spectral clustering
#' 
#' @param D Square matrix or `dist` object containing Euclidean distances 
#' between data points.
#' @param nn Integer specifying number of nearest neighbours to consider
#' @param nEig Integer specifying number of eigenvectors to retain.
#' @author Adapted by MRS from script by
#' [Nura Kawa](https://rpubs.com/nurakawa/spectral-clustering)
#' @examples
#' library('TreeTools', quietly = TRUE, warn.conflict = FALSE)
#' trees <- as.phylo(0:10, nTip = 8)
#' distances <- ClusteringInfoDistance(trees)
#' SpectralClustering(distances)
#' @family tree space functions
#' @export
SpectralClustering <- function(D, nn = 10, nEig = 2) {
  MutualKnnGraph <- function(D, nn = 10) {
    D <- as.matrix(D)
    
    # intialize the knn matrix
    knn_mat <- matrix(0, nrow = nrow(D), ncol = nrow(D))
    
    # find the 10 nearest neighbors for each point
    for (i in seq_len(nrow(D))) {
      neighbor_index <- order(D[i, ])[2:(nn + 1)]
      knn_mat[i, ][neighbor_index] <- 1
    }
    
    # Now we note that i,j are neighbors iff K[i,j] = 1 or K[j,i] = 1
    knn_mat <- knn_mat + t(knn_mat) # find mutual knn
    
    knn_mat[knn_mat == 2] <- 1
    
    # Return:
    knn_mat
  }
  
  GraphLaplacian <- function(W) {
    stopifnot(nrow(W) == ncol(W))
    
    g = colSums(W) # degrees of vertices
    n = nrow(W)
    
    # Return:
    D_half = diag(1 / sqrt(g))
    diag(n) - D_half %*% W %*% D_half
  }
  
  W <- MutualKnnGraph(D) # 1. matrix of similarities
  L <- GraphLaplacian(W) # 2. compute graph laplacian
  ei <- eigen(L, symmetric = TRUE) # 3. Compute the eigenvectors and values of L
  
  # Return the eigenvectors of the n_eig smallest eigenvalues:
  ei$vectors[, nrow(L) - rev(seq_len(nEig))]
}
