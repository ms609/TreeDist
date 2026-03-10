
#' Calculate distances in parallel
#' 
#' Accelerate distance calculation by employing multiple \acronym{CPU} workers.
#' 
#' ## OpenMP (recommended for `ClusteringInfoDistance` / `MutualClusteringInfo`)
#' 
#' When the package is built with \acronym{OpenMP} support (the default on
#' Linux and Windows; optional on macOS), pairwise
#' [`ClusteringInfoDistance()`] / [`MutualClusteringInfo()`] calculations use
#' an efficient multi-threaded code path automatically — no cluster setup is
#' required.
#' 
#' The number of \acronym{OpenMP} threads is controlled by the standard
#' `"mc.cores"` option:
#' 
#' ```r
#' options(mc.cores = parallel::detectCores())  # use all available cores
#' options(mc.cores = 4L)                        # or a fixed number
#' ```
#' 
#' The default is `1` (single-threaded). The \acronym{OpenMP} path is
#' substantially faster than the R-cluster path for typical analysis sizes and
#' is preferred when available.
#' 
#' ## R parallel cluster (other metrics)
#' 
#' For metrics that do not yet have an \acronym{OpenMP} batch implementation
#' (e.g. [`RobinsonFoulds()`], [`MatchingSplitDistance()`]), "TreeDist"
#' parallelizes via [`parCapply()`] using a cluster stored in
#' `options("TreeDist-cluster")`.
#' 
#' `StartParallel()` calls `parallel::makeCluster()` and registers the cluster.
#' 
#' `SetParallel()` registers a pre-existing cluster.
#' 
#' `StopParallel()` stops the current TreeDist cluster and releases resources.
#'
#' Note that R-cluster parallelism carries a serialisation overhead of ~2–3 s,
#' so it is only beneficial for large problems (roughly > 500 trees at 50 tips).
#' For [`ClusteringInfoDistance()`] the \acronym{OpenMP} path is faster at
#' every problem size and is used automatically when no cluster is registered.
#' 
#' @param \dots Parameters to pass to [`makeCluster()`].
#' @param cl An existing cluster.
#' 
#' @examples
#' # OpenMP parallelism: just set mc.cores before calling distance functions.
#' options(mc.cores = 2L)
#' # ClusteringInfoDistance(trees)  # now uses 2 OpenMP threads
#' options(mc.cores = NULL)  # restore default (single-threaded)
#' 
#' if (interactive()) { # R cluster: only worthwhile for non-OpenMP metrics
#'   library("TreeTools", quietly = TRUE)
#'   nCores <- ceiling(parallel::detectCores() / 2)
#'   StartParallel(nCores) # Takes a few seconds to set up processes
#'   GetParallel()
#'   RobinsonFoulds(as.phylo(0:6, 100))
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
