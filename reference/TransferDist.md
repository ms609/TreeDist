# Transfer dissimilarity between phylogenetic trees

Compute the transfer dissimilarity between phylogenetic trees, as
defined by Takazawa et al. (2026) . The transfer dissimilarity uses the
transfer distance (Lemoine et al. 2018) to compare bipartitions,
providing a finer-grained measure than the Robinson–Foulds distance.
Each branch in each tree is scored by how many taxa must be moved to
match its closest counterpart in the other tree, and these scores are
summed.

## Usage

``` r
TransferDist(
  tree1,
  tree2 = NULL,
  scale = TRUE,
  normalize = FALSE,
  reportMatching = FALSE
)

TransferDistance(
  tree1,
  tree2 = NULL,
  scale = TRUE,
  normalize = FALSE,
  reportMatching = FALSE
)

TransferDistSplits(
  splits1,
  splits2,
  nTip = attr(splits1, "nTip"),
  scale = TRUE,
  reportMatching = FALSE
)
```

## Arguments

- tree1, tree2:

  Trees of class `phylo`, with leaves labelled identically, or lists of
  such trees to undergo pairwise comparison. Where implemented,
  `tree2 = NULL` will compute distances between each pair of trees in
  the list `tree1` using a fast algorithm based on Day (1985) .

- scale:

  Logical; if `TRUE` (default), use the scaled transfer dissimilarity.
  If `FALSE`, use the unscaled transfer dissimilarity.

- normalize:

  If a numeric value is provided, this will be used as a maximum value
  against which to rescale results. If `TRUE`, results will be rescaled
  against a maximum value calculated from the specified tree sizes and
  topology, as specified in the "Normalization" section below. If
  `FALSE`, results will not be rescaled.

- reportMatching:

  Logical specifying whether to return the clade matchings as an
  attribute of the score.

- splits1, splits2:

  Logical matrices where each row corresponds to a leaf, either listed
  in the same order or bearing identical names (in any sequence), and
  each column corresponds to a split, such that each leaf is identified
  as a member of the ingroup (`TRUE`) or outgroup (`FALSE`) of the
  respective split.

- nTip:

  (Optional) Integer specifying the number of leaves in each split.

## Value

`TransferDist()` returns an object of class `dist` (if `tree2` is
`NULL`), a numeric matrix (if both `tree1` and `tree2` are lists), or a
numeric value (for a single pair). If `reportMatching = TRUE`, the
return value carries `matching` and `pairScores` attributes.

## Details

The `scaled` variant divides each branch's contribution by its depth
minus one, giving equal weight to all branches regardless of their depth
(analogous to the Robinson–Foulds distance). The `unscaled` variant uses
raw transfer distances, giving more weight to deep branches. Neither
variant satisfies the triangle inequality for trees with six or more
tips.

## Normalization

When `normalize = TRUE`, the scaled transfer dissimilarity is divided by
`2 * (n - 3)`, placing it in the range \[0, 1\]. The unscaled version is
divided by the maximum possible unscaled dissimilarity (following
Takazawa et al. (2026) ).

## References

Day WHE (1985). “Optimal algorithms for comparing trees with labeled
leaves.” *Journal of Classification*, **2**(1), 7–28.
[doi:10.1007/BF01908061](https://doi.org/10.1007/BF01908061) .  
  
Lemoine F, Domelevo Entfellner J, Wilkinson E, Correia D, Dávila Felipe
M, De Oliveira T, Gascuel O (2018). “Renewing Felsenstein's phylogenetic
bootstrap in the era of big data.” *Nature*, **556**(7702), 452–456.
[doi:10.1038/s41586-018-0043-0](https://doi.org/10.1038/s41586-018-0043-0)
.  
  
Takazawa Y, Takeda A, Hayamizu M, Gascuel O (2026). “Outperforming the
majority-rule consensus tree using fine-grained dissimilarity measures.”
*bioRxiv*.
[doi:10.64898/2026.03.16.712085](https://doi.org/10.64898/2026.03.16.712085)
.

## See also

Other tree distances:
[`HierarchicalMutualInfo()`](https://ms609.github.io/TreeDist/reference/HierarchicalMutualInfo.md),
[`JaccardRobinsonFoulds()`](https://ms609.github.io/TreeDist/reference/JaccardRobinsonFoulds.md),
[`KendallColijn()`](https://ms609.github.io/TreeDist/reference/KendallColijn.md),
[`MASTSize()`](https://ms609.github.io/TreeDist/reference/MASTSize.md),
[`MatchingSplitDistance()`](https://ms609.github.io/TreeDist/reference/MatchingSplitDistance.md),
[`NNIDist()`](https://ms609.github.io/TreeDist/reference/NNIDist.md),
[`NyeSimilarity()`](https://ms609.github.io/TreeDist/reference/NyeSimilarity.md),
[`PathDist()`](https://ms609.github.io/TreeDist/reference/PathDist.md),
[`Robinson-Foulds`](https://ms609.github.io/TreeDist/reference/Robinson-Foulds.md),
[`SPRDist()`](https://ms609.github.io/TreeDist/reference/SPRDist.md),
[`TreeDistance()`](https://ms609.github.io/TreeDist/reference/TreeDistance.md)

## Examples

``` r
library(TreeTools)
TransferDist(BalancedTree(8), PectinateTree(8))
#> [1] 3
TransferDist(BalancedTree(8), PectinateTree(8), scale = FALSE)
#> [1] 4

# All-pairs
TransferDist(as.phylo(0:5, 8))
#>     1   2   3   4   5
#> 2 2.0                
#> 3 2.0 2.0            
#> 4 3.0 3.0 3.0        
#> 5 3.0 3.0 3.0 1.0    
#> 6 3.5 3.5 1.5 1.5 1.5
```
