# Hierarchical Mutual Information

Calculate the Hierarchical Mutual Information (HMI) between two trees,
following the recursive algorithm of Perotti et al. (2020) .

This function was written during a code sprint: its documentation and
test cases have not yet been carefully scrutinized, and its
implementation may change without notice. Please alert the maintainer to
any issues you encounter.

## Usage

``` r
HierarchicalMutualInfo(tree1, tree2 = NULL, normalize = FALSE)

HMI(tree1, tree2 = NULL, normalize = FALSE)

SelfHMI(tree)

EHMI(tree1, tree2, precision = 0.01, minResample = 36)

AHMI(tree1, tree2, Mean = max, precision = 0.01, minResample = 36)
```

## Arguments

- normalize:

  If `FALSE`, return the raw HMI, in bits. If `TRUE`, normalize to range
  \[0,1\] by dividing by `max(SelfHMI(tree1), SelfHMI(tree2))`. If a
  function, divide by `normalize(SelfHMI(tree1), SelfHMI(tree2))`.

- tree, tree1, tree2:

  An object that can be coerced to an
  [`HPart`](https://ms609.github.io/TreeDist/reference/HPart.md) object.

- precision:

  Numeric; Monte Carlo sampling will terminate once the relative
  standard error falls below this value.

- minResample:

  Integer specifying minimum number of Monte Carlo samples to conduct.
  Avoids early termination when sample size is too small to reliably
  estimate the standard error of the mean.

- Mean:

  Function by which to combine the self-information of the two input
  hierarchies, in order to normalize the HMI.

## Value

`HierarchicalMutualInfo()` returns a numeric value representing the
hierarchical mutual information between the input trees, in bits,
normalized as specified. Higher values indicate more shared hierarchical
structure.

`SelfHMI()` returns the hierarchical mutual information of a tree
compared with itself, i.e. its hierarchical entropy (HH).

`EHMI()` returns the expected HMI against a uniform shuffling of element
labels, estimated by performing Monte Carlo resampling on the same
hierarchical structure until the relative standard error of the estimate
falls below `precision`. The attributes of the returned object list the
variance (`var`), standard deviation (`sd`), standard error of the mean
(`sem`) and relative error (`relativeError`) of the estimate, and the
number of Monte Carlo samples used to obtain it (`samples`).

`AHMI()` returns the adjusted HMI, normalized such that zero corresponds
to the expected HMI given a random shuffling of elements on the same
hierarchical structure. The attribute `sem` gives the standard error of
the estimate.

## Details

`HierarchicalMutualInfo()` computes the hierarchical mutual content of
trees (Perotti et al. 2015; Perotti et al. 2020) , which accounts for
the non-independence of information represented by nested splits.

`tree` is converted to a set of hierarchical partitions, and the mutual
information (in bits) is computed recursively; the contribution of a
node is given by:

\$\$I(t,s) = \log_2(n\_{ts}) - \dfrac{H\_{us} + H\_{tv} -
H\_{uv}}{n\_{ts}} + \text{mean}(I\_{uv})\$\$

Where:

- \\n\_{ts}\\ is the number of common elements between partitions

- \\H\_{us}, H\_{tv}, H\_{uv}\\ are entropy terms from child comparisons

- \\I\_{uv}\\ is the recursive HMI for child pairs

`AHMI()` calculates the adjusted hierarchical mutual information:
\$\$\text{AHMI}(t, s) = \dfrac{I(t, s) - \hat{I}(t, s)}{
\text{mean}(H(t), H(s)) - \hat{I}(t, s)}\$\$ Where:

- \\I(t, s)\\ is the hierarchical mutual information between `tree1` and
  `tree2`

- \\\hat{I}(t, s)\\ is the expected HMI between `tree1` and `tree2`,
  estimated by Monte Carlo sampling

- \\H(t), H(s)\\ is the entropy (self-mutual information) of each tree

## References

Perotti JI, Almeira N, Saracco F (2020). “Towards a Generalization of
Information Theory for Hierarchical Partitions.” *Physical Review E*,
**101**(6), 062148.
[doi:10.1103/PhysRevE.101.062148](https://doi.org/10.1103/PhysRevE.101.062148)
.  
  
Perotti JI, Tessone CJ, Caldarelli G (2015). “Hierarchical Mutual
Information for the Comparison of Hierarchical Community Structures in
Complex Networks.” *Physical Review E - Statistical, Nonlinear, and Soft
Matter Physics*, **92**(6), 062825-1–062825-13.
[doi:10.1103/PhysRevE.92.062825](https://doi.org/10.1103/PhysRevE.92.062825)
.

## See also

Other tree distances:
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
library("TreeTools", quietly = TRUE)

tree1 <- BalancedTree(8)
tree2 <- PectinateTree(8)

# Calculate HMI between two trees
HierarchicalMutualInfo(tree1, tree2)
#> [1] 0.4822863

# HMI normalized against the mean information content of tree1 and tree2
HierarchicalMutualInfo(tree1, tree2, normalize = mean)
#> [1] 0.2411432

# Normalized HMI above is equivalent to:
HMI(tree1, tree2) / mean(SelfHMI(tree1), SelfHMI(tree2))
#> [1] 0.2411432
# Expected mutual info for this pair of hierarchies
EHMI(tree1, tree2, precision = 0.1)
#> [1] 0.3027541
#> attr(,"var")
#> [1] 0.002918266
#> attr(,"sd")
#> [1] 0.05402098
#> attr(,"sem")
#> [1] 0.009003497
#> attr(,"samples")
#> [1] 36
#> attr(,"relativeError")
#> [1] 0.04290379
# The adjusted HMI normalizes against this expectation
AHMI(tree1, tree2, precision = 0.1)
#> [1] 0.07400499
#> attr(,"sem")
#> [1] 0.004966383
```
