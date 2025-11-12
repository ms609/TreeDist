# Kendall–Colijn distance

Calculate the Kendall–Colijn tree distance, a measure related to the
path difference.

## Usage

``` r
KendallColijn(tree1, tree2 = NULL, Vector = KCVector)

KCVector(tree)

PathVector(tree)

SplitVector(tree)

KCDiameter(tree)
```

## Arguments

- tree1, tree2:

  Trees of class `phylo`, with leaves labelled identically, or lists of
  such trees to undergo pairwise comparison. Where implemented,
  `tree2 = NULL` will compute distances between each pair of trees in
  the list `tree1` using a fast algorithm based on Day (1985) .

- Vector:

  Function converting a tree to a numeric vector.

  `KCVector`, the default, returns the number of edges between the
  common ancestor of each pair of leaves and the root of the tree (per
  Kendall and Colijn 2016) .

  `PathVector` returns the number of edges between each pair of leaves
  (per Steel and Penny 1993) .

  `SplitVector` returns the size of the smallest split that contains
  each pair of leaves (per Smith 2022 ).

- tree:

  A tree of class [`phylo`](https://rdrr.io/pkg/ape/man/read.tree.html).

## Value

`KendallColijn()` returns an array of numerics providing the distances
between each pair of trees in `tree1` and `tree2`, or `splits1` and
`splits2`.

`KCDiameter()` returns the value of the Kendall & Colijn's (2016) metric
distance between two pectinate trees with *n* leaves ordered in the
opposite direction, which I suggest (without any attempt at a proof) may
be a useful proxy for the diameter (i.e. maximum value) of the K–C
metric.

## Details

The Kendall–Colijn distance works by measuring, for each pair of leaves,
the distance from the most recent common ancestor of those leaves and
the root node. For a given tree, this produces a vector of values
recording the distance-from-the-root of each most recent common ancestor
of each pair of leaves.

Two trees are compared by taking the Euclidean distance between the
respective vectors. This is calculated by taking the square root of the
sum of the squares of the differences between the vectors.

An analogous distance can be created from any vector representation of a
tree. The split size vector metric (Smith 2022) is an attempt to mimic
the Kendall Colijn metric in situations where the position of the root
should not be afforded special significance; and the path distance
(Steel and Penny 1993) is a familiar alternative whose underlying vector
measures the distance of the last common ancestor of each pair of leaves
from the leaves themselves, i.e. the length of the path from one leaf to
another.

None of these vector-based methods performs as well as other tree
distances in measuring similarities in the relationships implied by a
pair of trees (Smith 2020) ; in particular, the Kendall Colijn metric is
strongly influenced by tree balance, and may not be appropriate for a
suite of common applications (Smith 2022) .

## Functions

- `KCVector()`: Creates a vector that characterises a rooted tree, as
  described in Kendall and Colijn (2016) .

- `PathVector()`: Creates a vector reporting the number of edges between
  each pair of leaves, per the path metric of Steel and Penny (1993) .

- `SplitVector()`: Creates a vector reporting the smallest split
  containing each pair of leaves, per the metric proposed in
  Smith (2022) .

## References

Day WHE (1985). “Optimal algorithms for comparing trees with labeled
leaves.” *Journal of Classification*, **2**(1), 7–28.
[doi:10.1007/BF01908061](https://doi.org/10.1007/BF01908061) .  
  
Kendall M, Colijn C (2016). “Mapping phylogenetic trees to reveal
distinct patterns of evolution.” *Molecular Biology and Evolution*,
**33**(10), 2735–2743.
[doi:10.1093/molbev/msw124](https://doi.org/10.1093/molbev/msw124) .  
  
Smith MR (2020). “Information theoretic Generalized Robinson-Foulds
metrics for comparing phylogenetic trees.” *Bioinformatics*, **36**(20),
5007–5013.
[doi:10.1093/bioinformatics/btaa614](https://doi.org/10.1093/bioinformatics/btaa614)
.  
  
Smith MR (2022). “Robust analysis of phylogenetic tree space.”
*Systematic Biology*, **71**(5), 1255–1270.
[doi:10.1093/sysbio/syab100](https://doi.org/10.1093/sysbio/syab100) .  
  
Steel MA, Penny D (1993). “Distributions of tree comparison metrics—some
new results.” *Systematic Biology*, **42**(2), 126–141.
[doi:10.1093/sysbio/42.2.126](https://doi.org/10.1093/sysbio/42.2.126) .

## See also

[`treespace::treeDist`](https://CRAN.R-project.org/package=treespace/vignettes/introduction.html)
is a more sophisticated, if more cumbersome, implementation that
supports lambda \> 0, i.e. use of edge lengths in tree comparison.

Other tree distances:
[`HierarchicalMutualInfo()`](https://ms609.github.io/TreeDist/dev/reference/HierarchicalMutualInfo.md),
[`JaccardRobinsonFoulds()`](https://ms609.github.io/TreeDist/dev/reference/JaccardRobinsonFoulds.md),
[`MASTSize()`](https://ms609.github.io/TreeDist/dev/reference/MASTSize.md),
[`MatchingSplitDistance()`](https://ms609.github.io/TreeDist/dev/reference/MatchingSplitDistance.md),
[`NNIDist()`](https://ms609.github.io/TreeDist/dev/reference/NNIDist.md),
[`NyeSimilarity()`](https://ms609.github.io/TreeDist/dev/reference/NyeSimilarity.md),
[`PathDist()`](https://ms609.github.io/TreeDist/dev/reference/PathDist.md),
[`Robinson-Foulds`](https://ms609.github.io/TreeDist/dev/reference/Robinson-Foulds.md),
[`SPRDist()`](https://ms609.github.io/TreeDist/dev/reference/SPRDist.md),
[`TreeDistance()`](https://ms609.github.io/TreeDist/dev/reference/TreeDistance.md)

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
KendallColijn(TreeTools::BalancedTree(8), TreeTools::PectinateTree(8))
#> [1] 11.48913

set.seed(0)
KendallColijn(TreeTools::BalancedTree(8), lapply(rep(8, 3), ape::rtree))
#> [1] 9.591663 5.567764 9.949874
KendallColijn(lapply(rep(8, 4), ape::rtree))
#>          1        2        3
#> 2 7.280110                  
#> 3 7.874008 8.185353         
#> 4 4.795832 7.071068 7.681146

KendallColijn(lapply(rep(8, 4), ape::rtree), Vector = SplitVector)
#>           1         2         3
#> 2 10.862780                    
#> 3 10.148892 12.124356          
#> 4  8.000000 11.489125  8.185353

# Notice that changing tree shape close to the root results in much
# larger differences
tree1 <- ape::read.tree(text = "(a, (b, (c, (d, (e, (f, (g, h)))))));")
tree2 <- ape::read.tree(text = "(a, ((b, c), (d, (e, (f, (g, h))))));")
tree3 <- ape::read.tree(text = "(a, (b, (c, (d, (e, ((f, g), h))))));")
trees <- c(tree1, tree2, tree3)
KendallColijn(trees)
#>          1        2
#> 2 4.000000         
#> 3 1.414214 4.242641
KendallColijn(trees, Vector = SplitVector)
#>          1        2
#> 2 2.449490         
#> 3 2.449490 3.162278
KCDiameter(4)
#> [1] 3.162278
KCDiameter(trees)
#> [1] 15.87451 15.87451 15.87451
```
