#### #'  Prepare a cluster for use in tree search functions
#### #' 
#### #' \code{PrepareCluster} creates a cluster of multiple cores and prepares it to analyse phylogenetic trees.
#### #'
#### #' @usage PrepareCluster(cores)
#### #' 
#### #' @param cores Number of cores to include in the cluster. 
#### #' 
#### #' @return Returns a reference to a cluster.
#### #' 
#### #' @author Martin R. Smith
#### #' 
#### #' @examples
#### #' \dontrun{
#### #'   PrepareCluster(4)
#### #' }
#### #' 
#### #' @importFrom parallel setDefaultCluster clusterEvalQ makeCluster
#### #' @export
#### PrepareCluster <- function (cores) {
####   cl <- makeCluster(getOption("cl.cores", cores))
####   clusterEvalQ(cl, {library(TreeSearch); NULL})
####   setDefaultCluster(cl)
####   attr(cl, 'cores') <- cores;
####   return(cl)
#### }