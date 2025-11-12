# Faithfulness of mapped distances

`MappingQuality()` calculates the trustworthiness and continuity of
mapped distances (Venna and Kaski 2001; Kaski et al. 2003) .
Trustworthiness measures, on a scale from 0–1, the degree to which
points that are nearby in a mapping are truly close neighbours;
continuity, the extent to which points that are truly nearby retain
their close spatial proximity in a mapping.

## Usage

``` r
MappingQuality(original, mapped, neighbours = 10L)

ProjectionQuality(original, mapped, neighbours = 10L)
```

## Arguments

- original, mapped:

  Square matrix or `dist` object containing original / mapped pairwise
  distances.

- neighbours:

  Integer specifying number of nearest neighbours to use in calculation.
  This should typically be small relative to the number of points.

## Value

`MappingQuality()` returns a named vector of length four, containing the
entries: `Trustworthiness`, `Continuity`, `TxC` (the product of these
values), and `sqrtTxC` (its square root).

## References

Kaski S, Nikkila J, Oja M, Venna J, Toronen P, Castren E (2003).
“Trustworthiness and metrics in visualizing similarity of gene
expression.” *BMC Bioinformatics*, **4**, 48.
[doi:10.1186/1471-2105-4-48](https://doi.org/10.1186/1471-2105-4-48) .  
  
Venna J, Kaski S (2001). “Neighborhood preservation in nonlinear
projection methods: an experimental study.” In Dorffner G, Bischof H,
Hornik K (eds.), *Artificial Neural Networks — ICANN 2001*, Lecture
Notes in Computer Science, 485–491.
[doi:10.1007/3-540-44668-0_68](https://doi.org/10.1007/3-540-44668-0_68)
.

## See also

Other tree space functions:
[`Islands()`](https://ms609.github.io/TreeDist/dev/reference/Islands.md),
[`MSTSegments()`](https://ms609.github.io/TreeDist/dev/reference/MSTSegments.md),
[`MapTrees()`](https://ms609.github.io/TreeDist/dev/reference/MapTrees.md),
[`SpectralEigens()`](https://ms609.github.io/TreeDist/dev/reference/SpectralEigens.md),
[`cluster-statistics`](https://ms609.github.io/TreeDist/dev/reference/cluster-statistics.md),
[`median.multiPhylo()`](https://ms609.github.io/TreeDist/dev/reference/median.multiPhylo.md)

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
library("TreeTools", quietly = TRUE)
trees <- as.phylo(0:10, nTip = 8)
distances <- ClusteringInfoDistance(trees)
mapping <- cmdscale(distances)
MappingQuality(distances, dist(mapping), 4)
#> Trustworthiness      Continuity             TxC         sqrtTxC 
#>       0.7929293       0.8737374       0.6928120       0.8323533 
```
