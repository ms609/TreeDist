# Per-pair transfer dissimilarity

Per-pair transfer dissimilarity

## Usage

``` r
cpp_transfer_dist(x, y, nTip)
```

## Arguments

- x, y:

  Raw matrices representing splits (from as.Splits()).

- nTip:

  Integer: number of tips.

## Value

A list with components:

- score_scaled: scaled transfer dissimilarity (double)

- score_unscaled: unscaled transfer dissimilarity (double)

- `matching_xy`: integer vector, best match in y for each split in x
  (1-based, NA if sentinel)

- `matching_yx`: integer vector, best match in x for each split in y
  (1-based, NA if sentinel)
