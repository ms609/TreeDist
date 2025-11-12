# Path distance

Calculate the path distance between rooted or unrooted trees.

## Usage

``` r
PathDist(tree1, tree2 = NULL)
```

## Arguments

- tree1, tree2:

  Trees of class `phylo`, with leaves labelled identically, or lists of
  such trees to undergo pairwise comparison. Where implemented,
  `tree2 = NULL` will compute distances between each pair of trees in
  the list `tree1` using a fast algorithm based on Day (1985) .

## Value

`PathDist()` returns a vector or distance matrix of distances between
trees.

## Details

This function is a faster alternative to the function
[`path.dist()`](https://klausvigo.github.io/phangorn/reference/treedist.html)
in the phangorn package, which can crash if the internal representation
of trees does not conform to certain (unspecified) expectations, and
which treats all trees as unrooted.

The path distance is calculated by tabulating the cladistic difference
(= topological distance) between each pair of tips in each tree. A
precursor to the path distance (Farris 1969) took the mean squared
difference between the elements of each tree's tabulation (Farris,
1973); the method used here is that proposed by Steel and Penny (1993) ,
which takes the square root of this sum. Other precursor measures are
described in Williams and Clifford (1971) and Phipps (1971) .

If a root node is present, trees are treated as rooted. To avoid
counting the root edge twice, use `UnrootTree(tree)` before calculating
the path distance.

Use of the path distance is discouraged as it emphasizes shallow
relationships at the expense of deeper (and arguably more fundamental)
relationships (Farris 1973) .

## References

Day WHE (1985). “Optimal algorithms for comparing trees with labeled
leaves.” *Journal of Classification*, **2**(1), 7–28.
[doi:10.1007/BF01908061](https://doi.org/10.1007/BF01908061) .  
  
Farris JS (1969). “A successive approximations approach to character
weighting.” *Systematic Biology*, **18**(4), 374–385.
[doi:10.2307/2412182](https://doi.org/10.2307/2412182) .  
  
Farris JS (1973). “On comparing the shapes of taxonomic trees.”
*Systematic Zoology*, **22**(1), 50–54.
[doi:10.2307/2412378](https://doi.org/10.2307/2412378) .  
  
Phipps JB (1971). “Dendrogram topology.” *Systematic Zoology*,
**20**(3), 306. [doi:10.2307/2412343](https://doi.org/10.2307/2412343)
.  
  
Steel MA, Penny D (1993). “Distributions of tree comparison metrics—some
new results.” *Systematic Biology*, **42**(2), 126–141.
[doi:10.1093/sysbio/42.2.126](https://doi.org/10.1093/sysbio/42.2.126)
.  
  
Williams WT, Clifford HT (1971). “On the comparison of two
classifications of the same set of elements.” *Taxon*, **20**(4),
519–522. [doi:10.2307/1218253](https://doi.org/10.2307/1218253) .

## See also

Other tree distances:
[`HierarchicalMutualInfo()`](https://ms609.github.io/TreeDist/reference/HierarchicalMutualInfo.md),
[`JaccardRobinsonFoulds()`](https://ms609.github.io/TreeDist/reference/JaccardRobinsonFoulds.md),
[`KendallColijn()`](https://ms609.github.io/TreeDist/reference/KendallColijn.md),
[`MASTSize()`](https://ms609.github.io/TreeDist/reference/MASTSize.md),
[`MatchingSplitDistance()`](https://ms609.github.io/TreeDist/reference/MatchingSplitDistance.md),
[`NNIDist()`](https://ms609.github.io/TreeDist/reference/NNIDist.md),
[`NyeSimilarity()`](https://ms609.github.io/TreeDist/reference/NyeSimilarity.md),
[`Robinson-Foulds`](https://ms609.github.io/TreeDist/reference/Robinson-Foulds.md),
[`SPRDist()`](https://ms609.github.io/TreeDist/reference/SPRDist.md),
[`TreeDistance()`](https://ms609.github.io/TreeDist/reference/TreeDistance.md)

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
library("TreeTools")

# Treating the two edges to the root node as distinct
PathDist(BalancedTree(7), PectinateTree(7))
#> [1] 5.656854

# Counting those two edges once
PathDist(UnrootTree(BalancedTree(7)), UnrootTree(PectinateTree(7)))
#> [1] 4.690416

PathDist(BalancedTree(7), as.phylo(0:2, 7))
#> [1] 8.944272 8.124038 8.944272
PathDist(as.phylo(0:2, 7), PectinateTree(7))
#> [1] 7.745967 8.246211 8.124038

PathDist(list(bal = BalancedTree(7), pec = PectinateTree(7)),
        as.phylo(0:2, 7))
#>          [,1]     [,2]     [,3]
#> [1,] 8.944272 8.124038 8.944272
#> [2,] 7.745967 8.246211 8.124038

PathDist(as.phylo(30:33, 8))
#>          1        2        3
#> 2 3.872983                  
#> 3 3.872983 3.464102         
#> 4 5.567764 6.928203 6.928203
 
```
