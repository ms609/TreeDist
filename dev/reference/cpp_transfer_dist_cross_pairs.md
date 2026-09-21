# Cross-pairs transfer dissimilarity (OpenMP)

Cross-pairs transfer dissimilarity (OpenMP)

## Usage

``` r
cpp_transfer_dist_cross_pairs(splits_a, splits_b, n_tip, scale, n_threads)
```

## Arguments

- splits_a, splits_b:

  Lists of raw matrices.

- n_tip:

  Number of tips.

- scale:

  Logical: use scaled transfer dissimilarity?

- n_threads:

  Number of OpenMP threads.

## Value

Numeric matrix of dimension `nA` x `nB`.
