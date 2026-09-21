# Pairwise mutual clustering information — batch computation

Internal function. Computes all pairwise MCI scores for a set of trees,
using OpenMP threads when available (falling back to single-threaded
execution otherwise). No interrupt checking is performed inside the
parallel region; the outer R call remains interruptible between batches.

## Usage

``` r
cpp_mutual_clustering_all_pairs(splits_list, n_tip, n_threads = 1L)
```

## Arguments

- splits_list:

  A list of split matrices (class `Splits` or `RawMatrix`), one per
  tree, all covering the same tip set. Typically the object returned by
  `as.Splits(trees, tipLabels = labs, asSplits = FALSE)`.

- n_tip:

  Integer; number of tips shared by all trees.

## Value

Numeric vector of length `n*(n-1)/2` containing pairwise MCI scores in
`combn(n, 2)` column-major order (i.e. the data payload of an R `dist`
object).
