# Median of a set of trees

Calculate the single binary tree that represents the geometric median –
an "average" – of a forest of tree topologies.

## Usage

``` r
# S3 method for class 'multiPhylo'
median(
  x,
  na.rm = FALSE,
  Distance = ClusteringInfoDistance,
  index = FALSE,
  breakTies = TRUE,
  ...
)
```

## Arguments

- x:

  Object of class `multiPhylo` containing phylogenetic trees.

- na.rm, ...:

  Unused; included for consistency with default function..

- Distance:

  Function to calculate distances between each pair of trees in `x`.

- index:

  Logical: if `TRUE`, return the index of the median tree(s); if
  `FALSE`, return the tree itself.

- breakTies:

  Logical: if `TRUE`, return a single tree with the minimum score; if
  `FALSE`, return all tied trees.

## Value

[`median()`](https://rdrr.io/r/stats/median.html) returns an object of
class `phylo` corresponding to the geometric median of a set of trees:
that is, the tree whose average distance from all other trees in the set
is lowest. If multiple trees tie in their average distance, the first
will be returned, unless `breakTies = FALSE`, in which case an object of
class `multiPhylo` containing all such trees will be returned.

## Details

The geometric median is the tree that exhibits the shortest average
distance from each other tree topology in the set. It represents an
"average" of a set of trees, though note that an unsampled tree may be
closer to the geometric "centre of gravity" of the input set – such a
tree would not be considered.

The result will depend on the metric chosen to calculate distances
between tree topologies. In the absence of a natural metric of tree
topologies, the default choice is
[`ClusteringInfoDistance()`](https://ms609.github.io/TreeDist/reference/TreeDistance.md)
– which discards branch length information. If specifying a different
function, be sure that it returns a difference, rather than a
similarity.

## See also

Consensus methods:
[`ape::consensus()`](https://rdrr.io/pkg/ape/man/consensus.html),
[`TreeTools::ConsensusWithout()`](https://ms609.github.io/TreeTools/reference/ConsensusWithout.html)

Other tree space functions:
[`Islands()`](https://ms609.github.io/TreeDist/reference/Islands.md),
[`MSTSegments()`](https://ms609.github.io/TreeDist/reference/MSTSegments.md),
[`MapTrees()`](https://ms609.github.io/TreeDist/reference/MapTrees.md),
[`MappingQuality()`](https://ms609.github.io/TreeDist/reference/MappingQuality.md),
[`SpectralEigens()`](https://ms609.github.io/TreeDist/reference/SpectralEigens.md),
[`cluster-statistics`](https://ms609.github.io/TreeDist/reference/cluster-statistics.md)

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
library("TreeTools", quietly = TRUE)
tenTrees <- as.phylo(1:10, nTip = 8)

# Default settings:
median(tenTrees)
#> 
#> Phylogenetic tree with 8 tips and 7 internal nodes.
#> 
#> Tip labels:
#>   t1, t2, t3, t4, t5, t6, ...
#> 
#> Rooted; no branch length.

# Robinson-Foulds distances include ties:
median(tenTrees, Distance = RobinsonFoulds, breakTies = FALSE)
#> 4 phylogenetic trees

# Be sure to use a distance function, rather than a similarity:
NyeDistance <- function(...) NyeSimilarity(..., similarity = FALSE)
median(tenTrees, Distance = NyeDistance)
#> 
#> Phylogenetic tree with 8 tips and 7 internal nodes.
#> 
#> Tip labels:
#>   t1, t2, t3, t4, t5, t6, ...
#> 
#> Rooted; no branch length.

# To analyse a list of trees that is not of class multiPhylo:
treeList <- lapply(1:10, as.phylo, nTip = 8)
class(treeList)
#> [1] "list"
median(structure(treeList, class = "multiPhylo"))
#> 
#> Phylogenetic tree with 8 tips and 7 internal nodes.
#> 
#> Tip labels:
#>   t1, t2, t3, t4, t5, t6, ...
#> 
#> Rooted; no branch length.
```
