# Transfer consensus (C++ implementation)

Transfer consensus (C++ implementation)

## Usage

``` r
cpp_transfer_consensus(
  splits_list,
  n_tip,
  scale,
  greedy_best_flag,
  init_majority,
  n_threads = 1L
)
```

## Arguments

- splits_list:

  List of raw matrices (one per tree), each from as.Splits().

- n_tip:

  Number of tips.

- scale:

  Logical: use scaled transfer distance?

- greedy_best_flag:

  Logical: TRUE for "best", FALSE for "first".

- init_majority:

  Logical: TRUE to start from majority-rule splits.

## Value

A `LogicalVector` of length n_splits indicating which pooled splits are
included in the consensus, plus attributes "raw_splits" (a raw matrix of
all unique splits) and "light_side" (integer vector).
