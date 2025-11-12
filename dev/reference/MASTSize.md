# Maximum Agreement Subtree size

Calculate the size or phylogenetic information content (Steel and Penny
2006) of the maximum agreement subtree between two phylogenetic trees,
i.e. the largest tree that can be obtained from both `tree1` and `tree2`
by deleting, but not rearranging, leaves, using the algorithm of
Valiente (2009) .

## Usage

``` r
MASTSize(tree1, tree2 = tree1, rooted = TRUE)

MASTInfo(tree1, tree2 = tree1, rooted = TRUE)
```

## Arguments

- tree1, tree2:

  Trees of class `phylo`, or lists of such trees to undergo pairwise
  comparison.

- rooted:

  Logical specifying whether to treat the trees as rooted.

## Value

`MASTSize()` returns an integer specifying the number of leaves in the
maximum agreement subtree.

`MASTInfo()` returns a vector or matrix listing the phylogenetic
information content, in bits, of the maximum agreement subtree.

## Details

Implemented for trees with up to 4096 tips. Contact the maintainer if
you need to process larger trees.

## References

Steel MA, Penny D (2006). “Maximum parsimony and the phylogenetic
information in multistate characters.” In Albert VA (ed.), *Parsimony,
Phylogeny, and Genomics*, 163–178. Oxford University Press, Oxford.  
  
Valiente G (2009). *Combinatorial Pattern Matching Algorithms in
Computational Biology using Perl and R*, CRC Mathematical and Computing
Biology Series. CRC Press, Boca Raton.

## See also

[`phangorn::mast()`](https://klausvigo.github.io/phangorn/reference/mast.html),
a slower implementation that also lists the leaves contained within the
subtree.

Other tree distances:
[`HierarchicalMutualInfo()`](https://ms609.github.io/TreeDist/dev/reference/HierarchicalMutualInfo.md),
[`JaccardRobinsonFoulds()`](https://ms609.github.io/TreeDist/dev/reference/JaccardRobinsonFoulds.md),
[`KendallColijn()`](https://ms609.github.io/TreeDist/dev/reference/KendallColijn.md),
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
 # for as.phylo, BalancedTree, PectinateTree:
library("TreeTools", quietly = TRUE)

MASTSize(PectinateTree(8), BalancedTree(8))
#> [1] 4
MASTInfo(PectinateTree(8), BalancedTree(8))
#> [1] 3.906891

MASTSize(BalancedTree(7), as.phylo(0:3, 7))
#> [1] 3 3 3 3
MASTSize(as.phylo(0:3, 7), PectinateTree(7))
#> [1] 6 5 5 5

MASTInfo(BalancedTree(7), as.phylo(0:3, 7))
#> [1] 1.584963 1.584963 1.584963 1.584963
MASTInfo(as.phylo(0:3, 7), PectinateTree(7))
#> [1] 9.884171 6.714246 6.714246 6.714246

MASTSize(list(Bal = BalancedTree(7), Pec = PectinateTree(7)),
         as.phylo(0:3, 7))
#>     [,1] [,2] [,3] [,4]
#> Bal    3    3    3    3
#> Pec    6    5    5    5
MASTInfo(list(Bal = BalancedTree(7), Pec = PectinateTree(7)),
         as.phylo(0:3, 7))
#>         [,1]     [,2]     [,3]     [,4]
#> Bal 1.584963 1.584963 1.584963 1.584963
#> Pec 9.884171 6.714246 6.714246 6.714246

CompareAll(as.phylo(0:4, 8), MASTSize)
#>   1 2 3 4 5
#> 1   7 7 7 6
#> 2 7   7 6 7
#> 3 7 7   6 6
#> 4 7 6 6   7
#> 5 6 7 6 7  
CompareAll(as.phylo(0:4, 8), MASTInfo)
#>           1         2         3         4         5
#> 1           13.343602 13.343602 13.343602  9.884171
#> 2 13.343602           13.343602  9.884171 13.343602
#> 3 13.343602 13.343602            9.884171  9.884171
#> 4 13.343602  9.884171  9.884171           13.343602
#> 5  9.884171 13.343602  9.884171 13.343602          
```
