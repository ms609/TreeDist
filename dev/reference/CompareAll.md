# Distances between each pair of trees

Calculate the distance between each tree in a list, and each other tree
in the same list.

## Usage

``` r
CompareAll(x, Func, FUN.VALUE = Func(x[[1]], x[[1]], ...), ...)
```

## Arguments

- x:

  List of trees, in the format expected by `Func()`.

- Func:

  distance function returning distance between two trees, e.g.
  [`path.dist()`](https://klausvigo.github.io/phangorn/reference/treedist.html).

- FUN.VALUE:

  Format of output of `Func()`, to be passed to
  [`vapply()`](https://rdrr.io/r/base/lapply.html). If unspecified,
  calculated by running `Func(x[[1]], x[[1]])`.

- ...:

  Additional parameters to pass to `Func()`.

## Value

`CompareAll()` returns a distance matrix of class `dist` detailing the
distance between each pair of trees. Identical trees are assumed to have
zero distance.

## Details

`CompareAll()` is not limited to tree comparisons: `Func` can be any
symmetric function.

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
# Generate a list of trees to compare
library("TreeTools", quietly = TRUE)
trees <- list(bal1 = BalancedTree(1:8), 
              pec1 = PectinateTree(1:8),
              pec2 = PectinateTree(c(4:1, 5:8)))

# Compare each tree with each other tree
CompareAll(trees, NNIDist)
#> $lower
#>      bal1 pec1 pec2
#> bal1         2    2
#> pec1    2         2
#> pec2    2    2     
#> 
#> $best_upper
#>      bal1 pec1 pec2
#> bal1         2    2
#> pec1    2         3
#> pec2    2    3     
#> 
#> $tight_upper
#>      bal1 pec1 pec2
#> bal1         2    2
#> pec1    2         3
#> pec2    2    3     
#> 
#> $loose_upper
#>      bal1 pec1 pec2
#> bal1         2    2
#> pec1    2         5
#> pec2    2    5     
#> 
#> $fack_upper
#>      bal1 pec1 pec2
#> bal1         2    2
#> pec1    2         5
#> pec2    2    5     
#> 
#> $li_upper
#>      bal1 pec1 pec2
#> bal1        10   10
#> pec1   10         8
#> pec2   10    8     
#> 
  
# Providing FUN.VALUE yields a modest speed gain:
dist <- CompareAll(trees, NNIDist, FUN.VALUE = integer(7))
#> Error in vapply(seq_along(i), function(k) {    cli_progress_update(1, .envir = parent.frame(2))    Func(i[[k]], j[[k]], ...)}, FUN.VALUE): values must be length 7,
#>  but FUN(X[[1]]) result is length 6
  
# View distances as a matrix
as.matrix(dist$lower)
#> Error in dist$lower: object of type 'closure' is not subsettable
```
