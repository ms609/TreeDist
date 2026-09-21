# Entropy in bits

Calculate the entropy of a vector of probabilities, in bits.
Probabilities should sum to one. Probabilities equalling zero will be
ignored.

## Usage

``` r
Entropy(...)

Ntropy(...)
```

## Arguments

- ...:

  Series of numerics, or single numeric vector, specifying probabilities
  of outcomes (for `Entropy()`) or counts (for `Ntropy()`).

## Value

`Entropy()` and `Ntropy()` return the entropy of the specified
probabilities or counts, in bits.

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
Entropy(1/2, 0, 1/2) # = 1
#> [1] 1
Entropy(rep(1/4, 4)) # = 2
#> [1] 2
Ntropy(c(2, 2, 0, 2, 2)) # = 2
#> [1] 2
```
