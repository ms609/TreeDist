
#' Calculate distances in parallel
#' 
#' Accelerate distance calculation by employing multiple \acronym{CPU} workers.
#' 
#' "TreeDist" parallelizes the calculation of tree to tree distances via
#' the [`parCapply()`] function, using a user-defined cluster specified in
#' `options("TreeDist-cluster")`.
#' 
#' `StartParallel()` calls `parallel::makeCluster()` and tells "TreeDist" to
#' use the created cluster.
#' 
#' `SetParallel()` tells "TreeDist" to use a pre-existing or user-specified 
#' cluster.
#' 
#' `StopParallel()` stops the current TreeDist cluster.
#' 
#' @param \dots Parameters to pass to [`makeCluster()`].
#' @param cl An existing cluster.
#' 
#' @examples
#' if (interactive()) { # Only run in terminal
#'   library("TreeTools", quietly = TRUE)
#'   nCores <- ceiling(detectCores() / 2)
#'   StartParallel(nCores) # Takes a few seconds to set up processes
#'   GetParallel()
#'   ClusteringInfoDistance(as.phylo(0:6, 100))
#'   StopParallel() # Returns system resources
#' }
#' @template MRS
#' @importFrom parallel makeCluster
#' @importFrom cli cli_alert_success cli_alert_danger
#' @export
StartParallel <- function(...) {
  cl <- makeCluster(...)
  cli_alert_success("Started cluster")
  options("TreeDist-cluster" = cl)
}

#' @rdname StartParallel
#' @return `StartParallel()` and `SetParallel()` return the previous value of
#' `options("TreeDist-cluster")`.
#' @export
SetParallel <- function(cl) {
  options("TreeDist-cluster" = cl)
}

#' @rdname StartParallel
#' @importFrom cli cli_alert_info
#' @return `GetParallel()` returns the currently specified cluster.
#' @export
GetParallel <- function(cl) {
  ret <- getOption("TreeDist-cluster")
  if (is.null(ret)) {
    cli_alert_info("No cluster currently specified")
  }
  ret
}

#' @rdname StartParallel
#' @param quietly Logical; if `TRUE`, do not warn when no cluster was running.
#' @export
#' @importFrom parallel stopCluster
#' @return `StopParallel()` returns `TRUE` if a cluster was destroyed,
#' `FALSE` otherwise.
StopParallel <- function(quietly = FALSE) {
  cluster <- getOption("TreeDist-cluster")
  if (!is.null(cluster)) {
    stopCluster(cluster)
    cli_alert_success("Cluster destroyed")
    options("TreeDist-cluster" = NULL)
    TRUE
  } else {
    if (!quietly) {
      cli_alert_danger("No cluster running")
    }
    FALSE
  }
}
