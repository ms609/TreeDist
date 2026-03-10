
#' Calculate distances in parallel
#' 
#' Accelerate distance calculation by employing multiple \acronym{CPU} workers.
#' 
#' ## OpenMP (recommended for all split-based metrics)
#' 
#' When the package is built with \acronym{OpenMP} support (the default on
#' Linux and Windows; optional on macOS), all pairwise split-based distance
#' calculations use an efficient multi-threaded batch path automatically —
#' no cluster setup is required.  The affected functions are:
#' 
#' - [`ClusteringInfoDistance()`] / [`MutualClusteringInfo()`]
#' - [`SharedPhylogeneticInfo()`] / [`DifferentPhylogeneticInfo()`]
#' - [`MatchingSplitInfo()`] / [`MatchingSplitInfoDistance()`]
#' - [`MatchingSplitDistance()`]
#' - [`InfoRobinsonFoulds()`]
#' - [`NyeSimilarity()`]
#' - [`JaccardRobinsonFoulds()`]
#' 
#' The number of \acronym{OpenMP} threads is controlled by the standard
#' `"mc.cores"` option:
#' 
#' ```r
#' options(mc.cores = parallel::detectCores())  # use all available cores
#' options(mc.cores = 4L)                        # or a fixed number
#' ```
#' 
#' The default is `1` (single-threaded).
#' 
#' ## R parallel cluster
#' 
#' `StartParallel()` creates an R socket cluster (via [`makeCluster()`]) and
#' registers it for use by TreeDist.  `SetParallel()` registers a pre-existing
#' cluster.  `StopParallel()` stops the cluster and releases resources.
#' 
#' **When to use `StartParallel()`:** for metrics that do not have an
#' \acronym{OpenMP} batch path, namely tree-object-based distances such as
#' [`NNIDist()`] and [`MASTSize()`] / [`MASTInfo()`], or any function called
#' via [`CompareAll()`].  R-cluster parallelism carries a serialisation overhead
#' of ~2–3 s, so it is only beneficial for large problems.
#' 
#' **When _not_ to use `StartParallel()`:** for the split-based metrics listed
#' above.  Registering a cluster disables the \acronym{OpenMP} batch path for
#' those functions, replacing a thread-local C++ loop with inter-process
#' communication — which is slower at every problem size measured.  Call
#' `StopParallel()` before computing split-based distances if a cluster is
#' active.
#' 
#' @param \dots Parameters to pass to [`makeCluster()`].
#' @param cl An existing cluster.
#' 
#' @examples
#' # OpenMP parallelism: set mc.cores before calling any split-based metric.
#' options(mc.cores = 2L)
#' # MutualClusteringInfo(trees)  # uses 2 OpenMP threads automatically
#' options(mc.cores = NULL)  # restore default (single-threaded)
#' 
#' if (interactive()) {
#'   # R cluster: beneficial for NNIDist, MASTSize/MASTInfo, CompareAll(), etc.
#'   # Do NOT activate while computing split-based distances (MCI, SPI, MSI, …)
#'   # as it bypasses the faster OpenMP path.
#'   library("TreeTools", quietly = TRUE)
#'   nCores <- ceiling(parallel::detectCores() / 2)
#'   StartParallel(nCores) # Takes a few seconds to set up processes
#'   GetParallel()
#'   CompareAll(as.phylo(0:6, 100), NNIDist)
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
