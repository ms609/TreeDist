# Create function to Efficiently compute distances from one tree to others in a set

`DistanceOracle()` prepares a set of trees and returns a function that
efficiently computes the distances from any one of them to all of them.

## Usage

``` r
DistanceOracle(trees, Func = ClusteringInfoDistance, ..., maxBytes = 2^31)
```

## Arguments

- trees:

  List of trees of class `phylo` or `Splits`, or a `multiPhylo` object.

- Func:

  Function giving the distance between trees: one of
  [`ClusteringInfoDistance()`](https://ms609.github.io/TreeDist/dev/reference/TreeDistance.md),
  [`InfoRobinsonFoulds()`](https://ms609.github.io/TreeDist/dev/reference/Robinson-Foulds.md),
  [`RobinsonFoulds()`](https://ms609.github.io/TreeDist/dev/reference/Robinson-Foulds.md),
  [`MatchingSplitDistance()`](https://ms609.github.io/TreeDist/dev/reference/MatchingSplitDistance.md),
  [`NyeSimilarity()`](https://ms609.github.io/TreeDist/dev/reference/NyeSimilarity.md)
  (returns a similarity unless `similarity = FALSE`), or
  [`JaccardRobinsonFoulds()`](https://ms609.github.io/TreeDist/dev/reference/JaccardRobinsonFoulds.md).

- ...:

  Further arguments to `Func`, such as `normalize`.

- maxBytes:

  Numeric giving the maximum memory, in bytes, that the prepared trees
  may occupy.

## Value

`DistanceOracle()` returns a function that takes a single argument `i`,
an integer vector that specifies the indices of one or more elements in
`trees`. That function returns the values of `Func` between `trees[[i]]`
and each tree in `trees`: a vector if `i` has length one, otherwise a
matrix with one column per entry of `i`.

## Details

In use cases that require one column of a distance matrix at a time, it
is faster to run a single preparation pass over all trees, than to
repeat the preparatory steps each time a new column is required.
`DistanceOracle()` produces a function that, when passed `i`, returns
the same as `Func(trees[[i]], trees, ...)`; repeated calls to that
function are faster than repeated calls to `Func`.

Columns are computed using `getOption("mc.cores")` threads.

## See also

Other tree distances:
[`HierarchicalMutualInfo()`](https://ms609.github.io/TreeDist/dev/reference/HierarchicalMutualInfo.md),
[`JaccardRobinsonFoulds()`](https://ms609.github.io/TreeDist/dev/reference/JaccardRobinsonFoulds.md),
[`KendallColijn()`](https://ms609.github.io/TreeDist/dev/reference/KendallColijn.md),
[`MASTSize()`](https://ms609.github.io/TreeDist/dev/reference/MASTSize.md),
[`MatchingSplitDistance()`](https://ms609.github.io/TreeDist/dev/reference/MatchingSplitDistance.md),
[`NNIDist()`](https://ms609.github.io/TreeDist/dev/reference/NNIDist.md),
[`NNIMaxStep()`](https://ms609.github.io/TreeDist/dev/reference/NNIMaxStep.md),
[`NyeSimilarity()`](https://ms609.github.io/TreeDist/dev/reference/NyeSimilarity.md),
[`PathDist()`](https://ms609.github.io/TreeDist/dev/reference/PathDist.md),
[`Robinson-Foulds`](https://ms609.github.io/TreeDist/dev/reference/Robinson-Foulds.md),
[`SPRDist()`](https://ms609.github.io/TreeDist/dev/reference/SPRDist.md),
[`TransferDist()`](https://ms609.github.io/TreeDist/dev/reference/TransferDist.md),
[`TreeDistance()`](https://ms609.github.io/TreeDist/dev/reference/TreeDistance.md)

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
library("TreeTools", quietly = TRUE)
trees <- as.phylo(0:9, nTip = 12)
Oracle <- DistanceOracle(trees, ClusteringInfoDistance)
Oracle(3)
#>  [1] 1.114993 1.114993 0.000000 2.417249 2.417249 1.375291 1.735396 2.329028
#>  [9] 1.995695 2.329028
all.equal(Oracle(3), ClusteringInfoDistance(trees[[3]], trees))
#> [1] TRUE
```
