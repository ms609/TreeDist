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

"TreeDist" parallelizes the calculation of tree to tree distances via
the `parCapply()` function, using a user-defined cluster specified in
`options("TreeDist-cluster")`.

`StartParallel()` calls
[`parallel::makeCluster()`](https://rdrr.io/r/parallel/makeCluster.html)
and tells "TreeDist" to use the created cluster.

`SetParallel()` tells "TreeDist" to use a pre-existing or user-specified
cluster.

`StopParallel()` stops the current TreeDist cluster.

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
if (interactive()) { # Only run in terminal
  library("TreeTools", quietly = TRUE)
  nCores <- ceiling(parallel::detectCores() / 2)
  StartParallel(nCores) # Takes a few seconds to set up processes
  GetParallel()
  ClusteringInfoDistance(as.phylo(0:6, 100))
  StopParallel() # Returns system resources
}
```
