# Generalized Robinson–Foulds distance

An internal function to calculate Generalized Robinson–Foulds distances
from splits.

## Usage

``` r
GeneralizedRF(
  splits1,
  splits2,
  nTip,
  PairScorer,
  maximize,
  reportMatching,
  ...
)
```

## Arguments

- splits1, splits2:

  Logical matrices where each row corresponds to a leaf, either listed
  in the same order or bearing identical names (in any sequence), and
  each column corresponds to a split, such that each leaf is identified
  as a member of the ingroup (`TRUE`) or outgroup (`FALSE`) of the
  respective split.

- nTip:

  Integer specifying the number of leaves in each split.

- PairScorer:

  function taking four arguments, `splits1`, `splits2`, `nSplits1`,
  `nSplits2`, which should return the score of each pair of splits in a
  two-dimensional matrix. Additional parameters may be specified via
  `...`.

- maximize:

  Logical specifying whether the optimal matching maximizes or minimizes
  the scores obtained by `PairScorer()`.

- reportMatching:

  Logical specifying whether to return the clade matchings as an
  attribute of the score.

- ...:

  Additional parameters to `PairScorer()`.

## Value

A numeric value specifying the score of the tree pairs under the
specified pair scorer. If `reportMatching = TRUE`, attribute also list:

- `matching`: which split in `splits2` is optimally matched to each
  split in `split1` (`NA` if not matched);

- `matchedSplits`: Textual representation of each match

- `matchedScores`: Scores for matched split.

- `pairScores`: Calculated scores for each possible matching of each
  split.

## Details

Note that no checks will be made to confirm that `splits1` and `splits2`
contain the same leaves in the same order. This is the responsibility of
the calling function.

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)
