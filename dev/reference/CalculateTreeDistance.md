# Wrapper for tree distance calculations

Calls tree distance functions from trees or lists of trees

## Usage

``` r
CalculateTreeDistance(Func, tree1, tree2 = NULL, reportMatching = FALSE, ...)
```

## Arguments

- Func:

  Tree distance function.

- tree1, tree2:

  Trees of class `phylo`, with leaves labelled identically, or lists of
  such trees to undergo pairwise comparison. Where implemented,
  `tree2 = NULL` will compute distances between each pair of trees in
  the list `tree1` using a fast algorithm based on Day (1985) .

- reportMatching:

  Logical specifying whether to return the clade matchings as an
  attribute of the score.

- ...:

  Additional arguments to `Func`.

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)
