# Approximate the Subtree Prune and Regraft (SPR) distance.

`SPRDist()` calculates an upper bound on the SPR distance between trees
using the heuristic method of de Oliveira Martins et al. (2008) . Other
approximations are available (e.g. Hickey et al. 2008, Goloboff 2008,
Whidden and Matsen 2018) .

## Usage

``` r
SPRDist(tree1, tree2 = NULL, method = "deOliveira", symmetric)

# S3 method for class 'phylo'
SPRDist(tree1, tree2 = NULL, method = "deOliveira", symmetric)

# S3 method for class 'list'
SPRDist(tree1, tree2 = NULL, method = "deOliveira", symmetric)

# S3 method for class 'multiPhylo'
SPRDist(tree1, tree2 = NULL, method = "deOliveira", symmetric)
```

## Arguments

- tree1, tree2:

  Trees of class `phylo`, with leaves labelled identically, or lists of
  such trees to undergo pairwise comparison. Where implemented,
  `tree2 = NULL` will compute distances between each pair of trees in
  the list `tree1` using a fast algorithm based on Day (1985) .

- method:

  Character specifying which method to use to approximate the SPR
  distance. Currently defaults to `"deOliveira"`, the only available
  option; a new method will eventually become the default.

- symmetric:

  Ignored (redundant after fix of
  [phangorn#97](https://github.com/KlausVigo/phangorn/issues/97)).

## Value

`SPRDist()` returns a vector or distance matrix of distances between
trees.

## References

Day WHE (1985). “Optimal algorithms for comparing trees with labeled
leaves.” *Journal of Classification*, **2**(1), 7–28.
[doi:10.1007/BF01908061](https://doi.org/10.1007/BF01908061) .  
  
de Oliveira Martins L, Leal E, Kishino H (2008). “Phylogenetic detection
of recombination with a Bayesian prior on the distance between trees.”
*PLoS One*, **3**(7), e2651.
[doi:10.1371/journal.pone.0002651](https://doi.org/10.1371/journal.pone.0002651)
.  
  
Goloboff PA (2008). “Calculating SPR distances between trees.”
*Cladistics*, **24**(4), 591-597.
[doi:10.1111/j.1096-0031.2007.00189.x](https://doi.org/10.1111/j.1096-0031.2007.00189.x)
.  
  
Hickey G, Dehne F, Rau-Chaplin A, Blouin C (2008). “SPR distance
computation for *unrooted* trees.” *Evolutionary Bioinformatics*, **4**,
EBO–S419. [doi:10.4137/EBO.S419](https://doi.org/10.4137/EBO.S419) .  
  
Whidden C, Matsen FA (2018). “Efficiently Inferring Pairwise Subtree
Prune-and-Regraft Adjacencies between Phylogenetic Trees.” *2018
Proceedings of the Meeting on Analytic Algorithmics and Combinatorics
(ANALCO)*, 77–91.
[doi:10.1137/1.9781611975062.8](https://doi.org/10.1137/1.9781611975062.8)
.

## See also

Exact calculation with
[TBRDist](https://ms609.github.io/TBRDist/reference/TreeRearrangementDistances.html)
functions `USPRDist()` and `ReplugDist()`.

phangorn function
[`SPR.dist()`](https://klausvigo.github.io/phangorn/reference/treedist.html)
employs the de Oliveira Martins et al. (2008) algorithm but can crash
when sent trees of certain formats, and tends to have a longer running
time.

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
[`TreeDistance()`](https://ms609.github.io/TreeDist/dev/reference/TreeDistance.md)

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
library("TreeTools", quietly = TRUE)

# Compare single pair of trees
SPRDist(BalancedTree(7), PectinateTree(7))
#> [1] 2

# Compare all pairs of trees        
SPRDist(as.phylo(30:33, 8))
#>   1 2 3 4
#> 1   1 1 1
#> 2 1   1 1
#> 3 1 1   1
#> 4 1 1 1  

# Compare each tree in one list with each tree in another
SPRDist(BalancedTree(7), as.phylo(0:2, 7))
#> [1] 2 2 2
SPRDist(as.phylo(0:2, 7), PectinateTree(7))
#> [1] 1 2 2

SPRDist(list(bal = BalancedTree(7), pec = PectinateTree(7)),
        as.phylo(0:2, 7))
#>     [,1] [,2] [,3]
#> bal    2    2    2
#> pec    1    2    2
```
