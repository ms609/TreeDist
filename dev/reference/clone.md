# Clone / duplicate an object `clone()` physically duplicates objects

Clone / duplicate an object `clone()` physically duplicates objects

## Usage

``` r
clone(x, ...)

# S3 method for class 'HPart'
clone(x, tipLabels = attr(x, "tip.label"), ...)
```

## Arguments

- x:

  the object to be cloned

- ...:

  additional parameters for methods

- tipLabels:

  Character vector specifying sequence in which to order tip labels.

## Value

`clone()` typically returns an object of the same class and "value" as
the input `x`.

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)
