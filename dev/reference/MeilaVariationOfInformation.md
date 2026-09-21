# Use variation of clustering information to compare pairs of splits

Compare a pair of splits viewed as clusterings of taxa, using the
variation of clustering information proposed by (Meila 2007) .

## Usage

``` r
MeilaVariationOfInformation(split1, split2)

MeilaMutualInformation(split1, split2)
```

## Arguments

- split1, split2:

  Logical vectors listing leaves in a consistent order, identifying each
  leaf as a member of the ingroup (`TRUE`) or outgroup (`FALSE`) of the
  split in question.

## Value

`MeilaVariationOfInformation()` returns the variation of (clustering)
information, measured in bits.

`MeilaMutualInformation()` returns the mutual information, measured in
bits.

## Details

This is equivalent to the mutual clustering information (Vinh et al.
2010) . For the total information content, multiply the VoI by the
number of leaves.

## References

Meila M (2007). “Comparing clusterings—an information based distance.”
*Journal of Multivariate Analysis*, **98**(5), 873–895.
[doi:10.1016/j.jmva.2006.11.013](https://doi.org/10.1016/j.jmva.2006.11.013)
.  
  
Vinh NX, Epps J, Bailey J (2010). “Information theoretic measures for
clusterings comparison: variants, properties, normalization and
correction for chance.” *Journal of Machine Learning Research*, **11**,
2837–2854.
[doi:10.1145/1553374.1553511](https://doi.org/10.1145/1553374.1553511) .

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
# Maximum variation = information content of each split separately
A <- TRUE
B <- FALSE
MeilaVariationOfInformation(c(A, A, A, B, B, B), c(A, A, A, A, A, A))
#> [1] 1
Entropy(c(3, 3) / 6) + Entropy(c(0, 6) / 6)
#> [1] 1

# Minimum variation = 0
MeilaVariationOfInformation(c(A, A, A, B, B, B), c(A, A, A, B, B, B))
#> [1] 0

# Not always possible for two evenly-sized splits to reach maximum
# variation of information
Entropy(c(3, 3) / 6) * 2  # = 2
#> [1] 2
MeilaVariationOfInformation(c(A, A, A,B ,B, B), c(A, B, A, B, A, B)) # < 2
#> [1] 1.836592

# Phylogenetically uninformative groupings contain spliting information
Entropy(c(1, 5) / 6)
#> [1] 0.6500224
MeilaVariationOfInformation(c(B, A, A, A, A, A), c(A, A, A, A, A, B))
#> [1] 1.203213
```
