# Largest clustering information distance reachable by one nearest neighbour interchange

`NNIMaxStep()` returns the largest [Clustering Information
Distance](https://ms609.github.io/TreeDist/dev/reference/TreeDistance.md)
that can separate an *n*-leaf tree from any tree that differs from it by
a single nearest neighbour interchange (NNI) move.

## Usage

``` r
NNIMaxStep(tree, normalize = FALSE)
```

## Arguments

- tree:

  Tree of class `phylo`, or list of trees of class `list` or
  `multiPhylo`, or an integer specifying the number of leaves in a tree.

- normalize:

  Logical specifying whether to normalize the distance against the
  summed clustering information of the two trees.

## Value

`NNIMaxStep()` returns a numeric vector, one entry per tree (or leaf
count), giving the largest attainable distance, in bits when
`normalize = FALSE`, or as a fraction in the range \[0, 1\] when
`normalize = TRUE`. `NA` is returned where *n* \< 4, as no NNI move
exists.

The vector bears attributes `"subtrees"`, giving the sizes of the four
subtrees around the moved edge, and `"splits"`, the sizes of the two
splits that the move exchanges. Where more than one tree is supplied,
each attribute is a list with one entry per tree.

## See also

The distance itself:
[`ClusteringInfoDistance()`](https://ms609.github.io/TreeDist/dev/reference/TreeDistance.md)

Diameter of the NNI metric:
[`NNIDiameter()`](https://ms609.github.io/TreeDist/dev/reference/NNIDist.md)

Other tree distances:
[`HierarchicalMutualInfo()`](https://ms609.github.io/TreeDist/dev/reference/HierarchicalMutualInfo.md),
[`JaccardRobinsonFoulds()`](https://ms609.github.io/TreeDist/dev/reference/JaccardRobinsonFoulds.md),
[`KendallColijn()`](https://ms609.github.io/TreeDist/dev/reference/KendallColijn.md),
[`MASTSize()`](https://ms609.github.io/TreeDist/dev/reference/MASTSize.md),
[`MatchingSplitDistance()`](https://ms609.github.io/TreeDist/dev/reference/MatchingSplitDistance.md),
[`NNIDist()`](https://ms609.github.io/TreeDist/dev/reference/NNIDist.md),
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
# Largest clustering information distance from a single NNI move
NNIMaxStep(8)  # exactly two bits for any multiple of four
#> [1] 2
#> attr(,"subtrees")
#> [1] 2 2 2 2
#> attr(,"splits")
#> [1] 4 4
NNIMaxStep(6)  # a little less otherwise
#> [1] 1.918296
#> attr(,"subtrees")
#> [1] 1 1 2 2
#> attr(,"splits")
#> [1] 2 3

# Read off the maximizing local topology
m6 <- NNIMaxStep(6)
attr(m6, "subtrees")
#> [1] 1 1 2 2
attr(m6, "splits")
#> [1] 2 3

# Vectorized over leaf counts
NNIMaxStep(4:8)
#> [1] 2.000000 1.901955 1.918296 1.929968 2.000000
#> attr(,"subtrees")
#> attr(,"subtrees")[[1]]
#> [1] 1 1 1 1
#> 
#> attr(,"subtrees")[[2]]
#> [1] 1 1 1 2
#> 
#> attr(,"subtrees")[[3]]
#> [1] 1 1 2 2
#> 
#> attr(,"subtrees")[[4]]
#> [1] 1 2 2 2
#> 
#> attr(,"subtrees")[[5]]
#> [1] 2 2 2 2
#> 
#> attr(,"splits")
#> attr(,"splits")[[1]]
#> [1] 2 2
#> 
#> attr(,"splits")[[2]]
#> [1] 2 3
#> 
#> attr(,"splits")[[3]]
#> [1] 2 3
#> 
#> attr(,"splits")[[4]]
#> [1] 3 3
#> 
#> attr(,"splits")[[5]]
#> [1] 4 4
#> 

# Computed for a given tree
library("TreeTools", quietly = TRUE)
NNIMaxStep(BalancedTree(19))
#> [1] 1.991546
#> attr(,"subtrees")
#> [1] 4 5 5 5
#> attr(,"splits")
#> [1] 9 9

# Normalized
NNIMaxStep(12, normalize = TRUE)
#> [1] 0.1460877
#> attr(,"subtrees")
#> [1] 3 3 3 3
#> attr(,"splits")
#> [1] 6 6
```
