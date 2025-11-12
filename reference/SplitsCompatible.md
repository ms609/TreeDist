# Are splits compatible?

Determine whether splits are compatible (concave); i.e. they can both
occur on a single tree.

## Usage

``` r
SplitsCompatible(split1, split2)
```

## Arguments

- split1, split2:

  Logical vectors listing leaves in a consistent order, identifying each
  leaf as a member of the ingroup (`TRUE`) or outgroup (`FALSE`) of the
  split in question.

## Value

`SplitsCompatible()` returns a logical specifying whether the splits
provided are compatible with one another.

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
A <- TRUE
B <- FALSE
SplitsCompatible(c(A, A, A, B, B, B),
                 c(A, A, B, B, B, B))
#> [1] TRUE
SplitsCompatible(c(A, A, A, B, B, B),
                 c(A, A, B, B, B, A))
#> [1] FALSE
```
