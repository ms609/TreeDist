# Calculate entropy of integer vector of counts

Wrapper for C++ function; no input checking is performed.
[`Ntropy()`](https://ms609.github.io/TreeDist/reference/Entropy.md) is
better suited for use where performance is not critical.

## Usage

``` r
entropy_int(n)
```

## Arguments

- n:

  a vector of integer counts

## Value

`entropy_int()` returns a numeric corresponding to the entropy of each
observation, in bits.
