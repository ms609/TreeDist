# Entropy of two splits

Calculate the entropy, joint entropy, entropy distance and information
content of two splits, treating each split as a division of *n* leaves
into two groups. Further details are available in a
[vignette](https://ms609.github.io/TreeDist/articles/information.html),
MacKay (2003) and Meila (2007) .

## Usage

``` r
SplitEntropy(split1, split2 = split1)
```

## Arguments

- split1, split2:

  Logical vectors listing leaves in a consistent order, identifying each
  leaf as a member of the ingroup (`TRUE`) or outgroup (`FALSE`) of the
  split in question.

## Value

A numeric vector listing, in bits:

- `H1` The entropy of split 1;

- `H2` The entropy of split 2;

- `H12` The joint entropy of both splits;

- `I` The mutual information of the splits;

- `Hd` The entropy distance (variation of information) of the splits.

## References

MacKay DJC (2003). *Information Theory, Inference, and Learning
Algorithms*. Cambridge University Press, Cambridge.
<https://www.inference.org.uk/itprnn/book.pdf>.  
  
Meila M (2007). “Comparing clusterings—an information based distance.”
*Journal of Multivariate Analysis*, **98**(5), 873–895.
[doi:10.1016/j.jmva.2006.11.013](https://doi.org/10.1016/j.jmva.2006.11.013)
.

## See also

Other information functions:
[`SplitSharedInformation()`](https://ms609.github.io/TreeDist/dev/reference/SplitSharedInformation.md),
[`TreeInfo`](https://ms609.github.io/TreeDist/dev/reference/TreeInfo.md)

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
A <- TRUE
B <- FALSE
SplitEntropy(c(A, A, A, B, B, B), c(A, A, B, B, B, B))
#>        H1        H2       H12         I        Hd 
#> 1.0000000 0.9182958 1.4591479 0.4591479 1.0000000 
```
