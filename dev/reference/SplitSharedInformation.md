# Shared information content of two splits

Calculate the phylogenetic information shared, or not shared, between
two splits. See the [accompanying
vignette](https://ms609.github.io/TreeDist/articles/information.html)
for definitions.

## Usage

``` r
SplitSharedInformation(n, A1, A2 = A1)

SplitDifferentInformation(n, A1, A2 = A1)

TreesConsistentWithTwoSplits(n, A1, A2 = A1)

LnTreesConsistentWithTwoSplits(n, A1, A2 = A1)

Log2TreesConsistentWithTwoSplits(n, A1, A2 = A1)

Log2TreesConsistentWithTwoSplits(n, A1, A2 = A1)
```

## Arguments

- n:

  Integer specifying the number of leaves

- A1, A2:

  Integers specifying the number of taxa in *A1* and *A2*, once the
  splits have been arranged such that *A1* fully overlaps with *A2*.

## Value

`TreesConsistentWithTwoSplits()` returns the number of unrooted
bifurcating trees consistent with two splits.

`SplitSharedInformation()` returns the phylogenetic information that two
splits have in common (Meila 2007) , in bits.

`SplitDifferentInformation()` returns the amount of phylogenetic
information distinct to one of the two splits, in bits.

## Details

Split *S1* divides *n* leaves into two splits, *A1* and *B1*. Split *S2*
divides the same leaves into the splits *A2* and *B2*.

Splits must be named such that *A1* fully overlaps with *A2*: that is to
say, all taxa in *A1* are also in *A2*, or *vice versa*. Thus, all taxa
in the smaller of *A1* and *A2* also occur in the larger.

## Functions

- `SplitDifferentInformation()`: Different information between two
  splits.

- `TreesConsistentWithTwoSplits()`: Number of trees consistent with two
  splits.

- `LnTreesConsistentWithTwoSplits()`: Natural logarithm of
  `TreesConsistentWithTwoSplits()`.

- `Log2TreesConsistentWithTwoSplits()`: Base two logarithm of
  `TreesConsistentWithTwoSplits()`.

- `Log2TreesConsistentWithTwoSplits()`: Base 2 logarithm of
  `TreesConsistentWithTwoSplits()`.

## References

Meila M (2007). “Comparing clusterings—an information based distance.”
*Journal of Multivariate Analysis*, **98**(5), 873–895.
[doi:10.1016/j.jmva.2006.11.013](https://doi.org/10.1016/j.jmva.2006.11.013)
.

## See also

Other information functions:
[`SplitEntropy()`](https://ms609.github.io/TreeDist/dev/reference/SplitEntropy.md),
[`TreeInfo`](https://ms609.github.io/TreeDist/dev/reference/TreeInfo.md)

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
  # Eight leaves, labelled A to H.
  # Split 1: ABCD|EFGH
  # Split 2: ABC|DEFGH
  # Let A1 = ABCD (four taxa), and A2 = ABC (three taxa).
  # A1 and A2 overlap (both contain ABC).
  
  TreesConsistentWithTwoSplits(n = 8, A1 = 4, A2 = 3)
#> [1] 45
  SplitSharedInformation(n = 8, A1 = 4, A2 = 3)
#> [1] 2.722466
  SplitDifferentInformation(n = 8, A1 = 4, A2 = 3)
#> [1] 5.129283

  # If splits are identical, then their shared information is the same
  # as the information of either split:
  SplitSharedInformation(n = 8, A1 = 3, A2 = 3)
#> [1] 5.044394
  TreeTools::SplitInformation(3, 5)
#> [1] 5.044394
```
