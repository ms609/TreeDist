# Calculate distances in parallel

Accelerate distance calculation by employing multiple CPU workers.

## Usage

``` r
StartParallel(...)

SetParallel(cl)

GetParallel(cl)

StopParallel(quietly = FALSE)
```

## Arguments

- ...:

  Parameters to pass to `makeCluster()`.

- cl:

  An existing cluster.

- quietly:

  Logical; if `TRUE`, do not warn when no cluster was running.

## Value

`StartParallel()` and `SetParallel()` return the previous value of
`options("TreeDist-cluster")`.

`GetParallel()` returns the currently specified cluster.

`StopParallel()` returns `TRUE` if a cluster was destroyed, `FALSE`
otherwise.

## Details

### OpenMP (recommended for all split-based metrics)

When the package is built with OpenMP support (the default on Linux and
Windows; optional on macOS), all pairwise split-based distance
calculations use an efficient multi-threaded batch path automatically —
no cluster setup is required. The affected functions are:

- [`ClusteringInfoDistance()`](https://ms609.github.io/TreeDist/dev/reference/TreeDistance.md)
  /
  [`MutualClusteringInfo()`](https://ms609.github.io/TreeDist/dev/reference/TreeDistance.md)

- [`SharedPhylogeneticInfo()`](https://ms609.github.io/TreeDist/dev/reference/TreeDistance.md)
  /
  [`DifferentPhylogeneticInfo()`](https://ms609.github.io/TreeDist/dev/reference/TreeDistance.md)

- [`MatchingSplitInfo()`](https://ms609.github.io/TreeDist/dev/reference/TreeDistance.md)
  /
  [`MatchingSplitInfoDistance()`](https://ms609.github.io/TreeDist/dev/reference/TreeDistance.md)

- [`MatchingSplitDistance()`](https://ms609.github.io/TreeDist/dev/reference/MatchingSplitDistance.md)

- [`InfoRobinsonFoulds()`](https://ms609.github.io/TreeDist/dev/reference/Robinson-Foulds.md)

- [`NyeSimilarity()`](https://ms609.github.io/TreeDist/dev/reference/NyeSimilarity.md)

- [`JaccardRobinsonFoulds()`](https://ms609.github.io/TreeDist/dev/reference/JaccardRobinsonFoulds.md)

The number of OpenMP threads is controlled by the standard `"mc.cores"`
option:

    options(mc.cores = parallel::detectCores())  # use all available cores
    options(mc.cores = 4L)                        # or a fixed number

The default is `1` (single-threaded).

### R parallel cluster

`StartParallel()` creates an R socket cluster (via `makeCluster()`) and
registers it for use by TreeDist. `SetParallel()` registers a
pre-existing cluster. `StopParallel()` stops the cluster and releases
resources.

**When to use `StartParallel()`:** for metrics that do not have an
OpenMP batch path, namely tree-object-based distances such as
[`NNIDist()`](https://ms609.github.io/TreeDist/dev/reference/NNIDist.md)
and
[`MASTSize()`](https://ms609.github.io/TreeDist/dev/reference/MASTSize.md)
/
[`MASTInfo()`](https://ms609.github.io/TreeDist/dev/reference/MASTSize.md),
or any function called via
[`CompareAll()`](https://ms609.github.io/TreeDist/dev/reference/CompareAll.md).
R-cluster parallelism carries a serialisation overhead of ~2–3 s, so it
is only beneficial for large problems.

**When *not* to use `StartParallel()`:** for the split-based metrics
listed above. Registering a cluster disables the OpenMP batch path for
those functions, replacing a thread-local C++ loop with inter-process
communication — which is slower at every problem size measured. Call
`StopParallel()` before computing split-based distances if a cluster is
active.

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
# OpenMP parallelism: set mc.cores before calling any split-based metric.
options(mc.cores = 2L)
# MutualClusteringInfo(trees)  # uses 2 OpenMP threads automatically
options(mc.cores = NULL)  # restore default (single-threaded)

if (interactive()) {
  # R cluster: beneficial for NNIDist, MASTSize/MASTInfo, CompareAll(), etc.
  # Do NOT activate while computing split-based distances (MCI, SPI, MSI, …)
  # as it bypasses the faster OpenMP path.
  library("TreeTools", quietly = TRUE)
  nCores <- ceiling(parallel::detectCores() / 2)
  StartParallel(nCores) # Takes a few seconds to set up processes
  GetParallel()
  CompareAll(as.phylo(0:6, 100), NNIDist)
  StopParallel() # Returns system resources
}
```
