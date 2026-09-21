# Normalize tree distances

`NormalizeInfo()` is an internal function used to normalize information
against a reference, such as the total information present in a pair of
trees.

## Usage

``` r
NormalizeInfo(
  unnormalized,
  tree1,
  tree2,
  InfoInTree,
  infoInBoth = NULL,
  how = TRUE,
  Combine = "+",
  ...
)
```

## Arguments

- unnormalized:

  Numeric value, vector or matrix to be normalized.

- tree1, tree2:

  Trees from which `unnormalized` was calculated.

- InfoInTree:

  Function to calculate the information content of each tree.

- infoInBoth:

  Optional numeric specifying information content of both trees
  independently. If unspecified (`NULL`), this will be calculated using
  the method specified by `how`.

- how:

  Method for normalization, perhaps specified using the `normalize`
  argument to a tree distance function. See details for options.

- ...:

  Additional parameters to `InfoInTree()` or `how()`.

## Value

`NormalizeInfo()` returns an object corresponding to the normalized
values of `unnormalized`.

## Details

The unnormalized value(s) are normalized by dividing by a denominator
calculated based on the `how` parameter. Valid options include:

- `FALSE`:

  No normalization is performed; the unnormalized values are returned.

- `TRUE`:

  Unless `infoInBoth` is specified, the information in each tree is
  computed using `InfoInTree()`, and the two values combined using
  `Combine()`.

- A numeric value, vector or matrix:

  `how` is used as the denominator; the returned value is
  `unnormalized / how`.

- A function:

  Unless `infoInBoth` is specified, the information in each tree is
  computed using `InfoInTree()`, and the two values combined using
  `how`. `NormalizeInfo(how = Func)` is thus equivalent to
  `NormalizeInfo(how = TRUE, Combine = Func)`.

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
library("TreeTools", quietly = TRUE)
pair1 <- c(BalancedTree(9), StarTree(9))
pair2 <- c(BalancedTree(9), PectinateTree(9))

# We'll let the number of nodes define the total information in a tree
Nnode(pair1)
#> [1] 8 1
Nnode(pair2)
#> [1] 8 8

# Let's normalize a unit distance
rawDist <- cbind(c(1, 1), c(1, 1))

# With `Combine = "+"`, the maximum distance is the sum of
# the information in each tree
denominator <- outer(Nnode(pair1), Nnode(pair2), "+")

NormalizeInfo(rawDist, pair1, pair2, InfoInTree = ape::Nnode, Combine = "+")
#>           [,1]      [,2]
#> [1,] 0.0625000 0.0625000
#> [2,] 0.1111111 0.1111111
rawDist / denominator
#>           [,1]      [,2]
#> [1,] 0.0625000 0.0625000
#> [2,] 0.1111111 0.1111111


# A denominator can be specified manually using `how`:
NormalizeInfo(rawDist, pair1, pair2, InfoInTree = ape::Nnode, how = 16)
#>        [,1]   [,2]
#> [1,] 0.0625 0.0625
#> [2,] 0.0625 0.0625
rawDist / 16
#>        [,1]   [,2]
#> [1,] 0.0625 0.0625
#> [2,] 0.0625 0.0625


# `how` also allows the denominator to be computed from trees:
outer(Nnode(pair1), Nnode(pair2), pmin)
#>      [,1] [,2]
#> [1,]    8    8
#> [2,]    1    1
NormalizeInfo(rawDist, pair1, pair2, InfoInTree = ape::Nnode, how = pmin)
#>       [,1]  [,2]
#> [1,] 0.125 0.125
#> [2,] 1.000 1.000
rawDist / outer(Nnode(pair1), Nnode(pair2), pmin)
#>       [,1]  [,2]
#> [1,] 0.125 0.125
#> [2,] 1.000 1.000
```
