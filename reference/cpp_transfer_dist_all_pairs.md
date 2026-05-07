# All-pairs transfer dissimilarity (OpenMP)

All-pairs transfer dissimilarity (OpenMP)

## Usage

``` r
cpp_transfer_dist_all_pairs(splits_list, n_tip, scale, n_threads)
```

## Arguments

- splits_list:

  List of raw matrices (one per tree).

- n_tip:

  Number of tips.

- scale:

  Logical: use scaled transfer dissimilarity?

- n_threads:

  Number of OpenMP threads.

## Value

Numeric vector of length choose(N,2) in dist order.
