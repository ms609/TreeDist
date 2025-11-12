# Eigenvalues for spectral clustering

Spectral clustering emphasizes nearest neighbours when forming clusters;
it avoids some of the issues that arise from clustering around means /
medoids.

## Usage

``` r
SpectralEigens(D, nn = 10L, nEig = 2L)

SpectralClustering(D, nn = 10L, nEig = 2L)
```

## Arguments

- D:

  Square matrix or `dist` object containing Euclidean distances between
  data points.

- nn:

  Integer specifying number of nearest neighbours to consider

- nEig:

  Integer specifying number of eigenvectors to retain.

## Value

`SpectralEigens()` returns spectral eigenvalues that can then be
clustered using a method of choice.

## See also

Other tree space functions:
[`Islands()`](https://ms609.github.io/TreeDist/dev/reference/Islands.md),
[`MSTSegments()`](https://ms609.github.io/TreeDist/dev/reference/MSTSegments.md),
[`MapTrees()`](https://ms609.github.io/TreeDist/dev/reference/MapTrees.md),
[`MappingQuality()`](https://ms609.github.io/TreeDist/dev/reference/MappingQuality.md),
[`cluster-statistics`](https://ms609.github.io/TreeDist/dev/reference/cluster-statistics.md),
[`median.multiPhylo()`](https://ms609.github.io/TreeDist/dev/reference/median.multiPhylo.md)

## Author

Adapted by MRS from script by [Nura
Kawa](https://rpubs.com/nurakawa/spectral-clustering)

## Examples

``` r
library("TreeTools", quietly = TRUE)
trees <- as.phylo(0:18, nTip = 8)
distances <- ClusteringInfoDistance(trees)
eigens <- SpectralEigens(distances)
# Perform clustering:
clusts <- KMeansPP(dist(eigens), k = 3)
plot(eigens, pch = 15, col = clusts$cluster)

plot(cmdscale(distances), pch = 15, col = clusts$cluster)
```
