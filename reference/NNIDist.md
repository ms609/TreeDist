# Approximate Nearest Neighbour Interchange distance

Use the approach of Li et al. (1996) to approximate the Nearest
Neighbour Interchange distance (Robinson 1971) between phylogenetic
trees.

## Usage

``` r
NNIDist(tree1, tree2 = tree1)

NNIDiameter(tree)
```

## Arguments

- tree1, tree2:

  Single trees of class `phylo` to undergo comparison.

- tree:

  Object of supported class representing a tree or list of trees, or an
  integer specifying the number of leaves in a tree/trees.

## Value

`NNIDist()` returns, for each pair of trees, a named vector containing
three integers:

- `lower` is a lower bound on the NNI distance, and corresponds to the
  RF distance between the trees.

- `tight_upper` is an upper bound on the distance, based on calculated
  maximum diameters for trees with \< 13 leaves. *NA* is returned if
  trees are too different to employ this approach.

- `loose_upper` is a looser upper bound on the distance, using *n* log
  *n* + O(*n*).

`NNIDiameter()` returns a matrix specifying (bounds on) the diameter of
the NNI distance metric on the specified tree(s). Columns correspond to:

- `liMin`: \$\$n - 3\$\$, a lower bound on the diameter (Li *et al.*
  1996);

- `fackMin`: Lower bound on diameter following Fack *et al*. (2002),
  i.e. \$\$\log2{N!} / 4\$\$;

- `min`: The larger of `liMin` and `fackMin`;

- `exact`: The exact value of the diameter, where *n* \< 13;

- `liMax`: Upper bound on diameter following Li *et al.* (1996), i.e.
  \$\$n \log2{n} + \textrm{O}(n)\$\$;

- `fackMax`: Upper bound on diameter following Fack *et al*. (2002),
  i.e. (\$\$N - 2\$\$) ceiling(\$\$\log2{n}\$\$)

  - *N*;

- `max`: The smaller of `liMax` and `fackMax`;

where *n* is the number of leaves, and *N* the number of internal nodes,
i.e. \$\$n - 2\$\$.

## Details

In brief, this approximation algorithm works by identifying edges in one
tree that do not match edges in the second. Each of these edges must
undergo at least one NNI operation in order to reconcile the trees.
Edges that match in both trees need never undergo an NNI operation, and
divide each tree into smaller regions. By "cutting" matched edges into
two, a tree can be divided into a number of regions that solely comprise
unmatched edges.

These regions can be viewed as separate trees that need to be
reconciled. One way to reconcile these trees is to conduct a series of
NNI operations that reduce a tree to a pectinate (caterpillar) tree,
then to conduct an analogue of the mergesort algorithm. This takes at
most *n* log *n* + O(*n*) NNI operations, and provides a loose upper
bound on the NNI score. The maximum number of moves for an *n*-leaf tree
([OEIS A182136](https://oeis.org/A182136)) can be calculated exactly for
small trees (Fack et al. 2002) ; this provides a tighter upper bound,
but is unavailable for *n* \> 12. `NNIDiameter()` reports the limits on
this bound.

|           |     |     |     |     |     |     |     |     |     |     |     |     |     |
|-----------|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|
| Leaves:   | 1   | 2   | 3   | 4   | 5   | 6   | 7   | 8   | 9   | 10  | 11  | 12  | 13  |
| Diameter: | 0   | 0   | 0   | 1   | 3   | 5   | 7   | 10  | 12  | 15  | 18  | 21  | ?   |

## References

Fack V, Lievens S, Van der Jeugt J (2002). “On the diameter of the
rotation graph of binary coupling trees.” *Discrete Mathematics*,
**245**(1-3), 1–18.
[doi:10.1016/S0012-365X(01)00418-6](https://doi.org/10.1016/S0012-365X%2801%2900418-6)
.  
  
Li M, Tromp J, Zhang L (1996). “Some notes on the nearest neighbour
interchange distance.” In Goos G, Hartmanis J, Leeuwen J, Cai J, Wong CK
(eds.), *Computing and Combinatorics*, volume 1090, 343–351. Springer,
Berlin, Heidelberg. ISBN 978-3-540-61332-9 978-3-540-68461-9,
[doi:10.1007/3-540-61332-3_168](https://doi.org/10.1007/3-540-61332-3_168)
.  
  
Robinson DF (1971). “Comparison of labeled trees with valency three.”
*Journal of Combinatorial Theory, Series B*, **11**(2), 105–119.
[doi:10.1016/0095-8956(71)90020-7](https://doi.org/10.1016/0095-8956%2871%2990020-7)
.

## See also

Other tree distances:
[`HierarchicalMutualInfo()`](https://ms609.github.io/TreeDist/reference/HierarchicalMutualInfo.md),
[`JaccardRobinsonFoulds()`](https://ms609.github.io/TreeDist/reference/JaccardRobinsonFoulds.md),
[`KendallColijn()`](https://ms609.github.io/TreeDist/reference/KendallColijn.md),
[`MASTSize()`](https://ms609.github.io/TreeDist/reference/MASTSize.md),
[`MatchingSplitDistance()`](https://ms609.github.io/TreeDist/reference/MatchingSplitDistance.md),
[`NyeSimilarity()`](https://ms609.github.io/TreeDist/reference/NyeSimilarity.md),
[`PathDist()`](https://ms609.github.io/TreeDist/reference/PathDist.md),
[`Robinson-Foulds`](https://ms609.github.io/TreeDist/reference/Robinson-Foulds.md),
[`SPRDist()`](https://ms609.github.io/TreeDist/reference/SPRDist.md),
[`TreeDistance()`](https://ms609.github.io/TreeDist/reference/TreeDistance.md)

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
library("TreeTools", quietly = TRUE)

NNIDist(BalancedTree(7), PectinateTree(7))
#>       lower  best_lower tight_upper  best_upper loose_upper  fack_upper 
#>           2           2           2           2           4           4 
#>    li_upper 
#>          10 

NNIDist(BalancedTree(7), as.phylo(0:2, 7))
#>             [,1] [,2] [,3]
#> lower          4    3    4
#> best_lower     7    5    7
#> tight_upper    7    5    7
#> best_upper     7    5    7
#> loose_upper   14    8   14
#> fack_upper    14    8   14
#> li_upper      16   13   16
NNIDist(as.phylo(0:2, 7), PectinateTree(7))
#>             [,1] [,2] [,3]
#> lower          4    4    4
#> best_lower     7    7    7
#> tight_upper    7    7    7
#> best_upper     7    7    7
#> loose_upper   14   14   14
#> fack_upper    14   14   14
#> li_upper      16   16   16

NNIDist(list(bal = BalancedTree(7), pec = PectinateTree(7)),
        as.phylo(0:2, 7))
#> , , lower
#> 
#>     [,1] [,2] [,3]
#> bal    4    3    4
#> pec    4    4    4
#> 
#> , , best_lower
#> 
#>     [,1] [,2] [,3]
#> bal    7    5    7
#> pec    7    7    7
#> 
#> , , tight_upper
#> 
#>     [,1] [,2] [,3]
#> bal    7    5    7
#> pec    7    7    7
#> 
#> , , best_upper
#> 
#>     [,1] [,2] [,3]
#> bal    7    5    7
#> pec    7    7    7
#> 
#> , , loose_upper
#> 
#>     [,1] [,2] [,3]
#> bal   14    8   14
#> pec   14   14   14
#> 
#> , , fack_upper
#> 
#>     [,1] [,2] [,3]
#> bal   14    8   14
#> pec   14   14   14
#> 
#> , , li_upper
#> 
#>     [,1] [,2] [,3]
#> bal   16   13   16
#> pec   16   16   16
#> 

CompareAll(as.phylo(30:33, 8), NNIDist)
#> $lower
#>   1 2 3 4
#> 1   1 1 2
#> 2 1   1 2
#> 3 1 1   2
#> 4 2 2 2  
#> 
#> $best_lower
#>   1 2 3 4
#> 1   1 1 3
#> 2 1   1 3
#> 3 1 1   3
#> 4 3 3 3  
#> 
#> $tight_upper
#>   1 2 3 4
#> 1   1 1 3
#> 2 1   1 3
#> 3 1 1   3
#> 4 3 3 3  
#> 
#> $best_upper
#>   1 2 3 4
#> 1   1 1 3
#> 2 1   1 3
#> 3 1 1   3
#> 4 3 3 3  
#> 
#> $loose_upper
#>   1 2 3 4
#> 1   2 2 5
#> 2 2   2 5
#> 3 2 2   5
#> 4 5 5 5  
#> 
#> $fack_upper
#>   1 2 3 4
#> 1   2 2 5
#> 2 2   2 5
#> 3 2 2   5
#> 4 5 5 5  
#> 
#> $li_upper
#>   1 2 3 4
#> 1   5 5 8
#> 2 5   5 8
#> 3 5 5   8
#> 4 8 8 8  
#> 
```
