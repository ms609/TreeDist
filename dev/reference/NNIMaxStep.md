# Largest clustering information distance reachable by one NNI move

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

  Object of supported class representing a tree or list of trees, or an
  integer specifying the number of leaves in a tree/trees.

- normalize:

  Logical specifying whether to normalize the distance against the
  summed clustering information of the two trees, giving the largest
  *normalized* Clustering Information Distance attainable by one NNI
  move.

## Value

`NNIMaxStep()` returns a numeric vector, one entry per tree (or leaf
count), giving the largest attainable distance – in bits when
`normalize = FALSE`, or as a fraction in the range \[0, 1\] when
`normalize = TRUE`. Entries are `NA` where *n* \< 4, as no NNI move
exists. Two attributes record the maximizing configuration:
`"subtrees"`, the sizes of the four subtrees around the moved edge, and
`"splits"`, the sizes of the two splits that the move exchanges. Where
more than one tree is supplied, each attribute is a list with one entry
per tree.

## Details

A single NNI move around an internal edge exchanges two of the four
subtrees that meet at that edge, changing exactly one split. Because
[`ClusteringInfoDistance()`](https://ms609.github.io/TreeDist/dev/reference/TreeDistance.md)
(Smith 2020) scores trees by an optimal matching of their splits, every
*unchanged* split matches its counterpart perfectly and contributes
nothing to the distance. The distance of an NNI move therefore depends
only on the old split *S* and the new split *S'*, and hence only on the
leaf counts \\a, b, c, d\\ of the four subtrees \\(a + b + c + d = n)\\.
Writing \\H\\ for entropy in bits, it equals their entropy distance
(variation of [clustering
information](https://ms609.github.io/TreeDist/dev/reference/TreeInfo.md),
(Meila 2007) ): \$\$d(a, b, c, d) = 2 H(a, b, c, d) - H(a + b) - H(a +
c).\$\$ The maximum can thus be found deterministically, by optimizing
over the four subtree sizes, rather than by sampling random trees. At
the optimum the four subtrees are as equal as possible (each \\\lfloor
n/4 \rfloor\\ or \\\lceil n/4 \rceil\\), so the maximizing topology is a
*local* property of a single edge, realized by *any* tree that contains
such an edge – not by a globally pectinate or balanced tree. The value
is exactly two bits when four divides *n* (four equal subtrees, whose
splits are then balanced and independent), and slightly less otherwise.
The exchanged splits themselves need not be balanced: at \\n = 6\\ the
optimal subtrees are \\\\2, 2, 1, 1\\\\ but the move runs between a \\4
\| 2\\ split and a \\3 \| 3\\ split.

## The normalized maximum is also local

With `normalize = TRUE` the distance is divided by the combined
clustering entropy of the two trees, \\CE(T) + CE(T')\\, which *does*
depend on the wider topology. The largest normalized move can
nevertheless be found without searching tree space, as follows.

\\T\\ and \\T'\\ share every split except the one the move changes, so
\$\$CE(T) + CE(T') = 2 E + H(S) + H(S'),\$\$ where \\E\\ is the summed
entropy of the \\n - 4\\ shared splits and \\H(S), H(S')\\ – the
entropies of the exchanged splits – are fixed by \\a, b, c, d\\. Each
shared split is cut by an edge lying *inside* one of the four subtrees,
so \\E\\ decomposes into four independent terms, \$\$E = c(a) + c(b) +
c(c) + c(d),\$\$ where \\c(s)\\ is the clustering entropy contributed by
a subtree of \\s\\ leaves. The numerator \\d(a, b, c, d)\\ is
independent of the subtree *shapes*, so for fixed sizes the ratio is
largest when each \\c(s)\\ is as *small* as possible – and each subtree
is minimized independently. The minimum \\c(s)\\ satisfies the recursion
\$\$c(1) = 0, \qquad c(s) = H(s) + \min\_{1 \le i \le s/2}\\c(i) + c(s -
i)\\,\$\$ computed once by dynamic programming. The largest normalized
NNI move is then \$\$\max\_{a, b, c, d} \frac{d(a, b, c, d)}{2\[c(a) +
c(b) + c(c) + c(d)\] + H(S) + H(S')},\$\$ a purely local optimization
over the four subtree sizes. (The minimizing subtree shape is not in
general balanced; for \\s = 6\\, for instance, it splits \\2 + 4\\
rather than \\3 + 3\\.)

## References

Meila M (2007). “Comparing clusterings—an information based distance.”
*Journal of Multivariate Analysis*, **98**(5), 873–895.
[doi:10.1016/j.jmva.2006.11.013](https://doi.org/10.1016/j.jmva.2006.11.013)
.  
  
Smith MR (2020). “Information theoretic Generalized Robinson-Foulds
metrics for comparing phylogenetic trees.” *Bioinformatics*, **36**(20),
5007–5013.
[doi:10.1093/bioinformatics/btaa614](https://doi.org/10.1093/bioinformatics/btaa614)
.

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
NNIMaxStep(8)  # exactly two bits: eight is a multiple of four
#> [1] 2
#> attr(,"subtrees")
#> [1] 2 2 2 2
#> attr(,"splits")
#> [1] 4 4
NNIMaxStep(6)  # a little less
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

# Vectorized over leaf counts, and accepting a tree
NNIMaxStep(4:12)
#> [1] 2.000000 1.901955 1.918296 1.929968 2.000000 1.967723 1.970951 1.973591
#> [9] 2.000000
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
#> attr(,"subtrees")[[6]]
#> [1] 2 2 2 3
#> 
#> attr(,"subtrees")[[7]]
#> [1] 2 2 3 3
#> 
#> attr(,"subtrees")[[8]]
#> [1] 2 3 3 3
#> 
#> attr(,"subtrees")[[9]]
#> [1] 3 3 3 3
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
#> attr(,"splits")[[6]]
#> [1] 4 5
#> 
#> attr(,"splits")[[7]]
#> [1] 4 5
#> 
#> attr(,"splits")[[8]]
#> [1] 5 5
#> 
#> attr(,"splits")[[9]]
#> [1] 6 6
#> 
library("TreeTools", quietly = TRUE)
NNIMaxStep(BalancedTree(19))
#> [1] 1.991546
#> attr(,"subtrees")
#> [1] 4 5 5 5
#> attr(,"splits")
#> [1] 9 9

# Normalized: solved locally, without searching tree space
NNIMaxStep(12, normalize = TRUE)
#> [1] 0.1460877
#> attr(,"subtrees")
#> [1] 3 3 3 3
#> attr(,"splits")
#> [1] 6 6
```
