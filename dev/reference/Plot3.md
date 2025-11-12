# Pseudo-3D plotting

`Plot3()` displays three-dimensional data in two dimensions, reflecting
the third dimension with point scaling, overlap and fogging. Points with
a lower `z` value are smaller than, fainter than, and overlapped by
points with a higher value.

## Usage

``` r
Plot3(
  x,
  y = NULL,
  z = NULL,
  pch = par("pch"),
  col = par("col"),
  bg = NA,
  cex = 1,
  axes = TRUE,
  frame.plot = axes,
  plot.bg = NA,
  fog = 1/2,
  shrink = 1/2,
  add = FALSE,
  ...
)
```

## Arguments

- x, y, z:

  Coordinates of points to plot.

- bg, cex, col, pch, add, axes, frame.plot, ...:

  Parameters passed to
  [`plot.default()`](https://rdrr.io/r/graphics/plot.default.html).

- plot.bg:

  Colour with which to fill plot area, used as fog colour.

- fog:

  Numeric from zero (no fading) to one (furthest points are invisible)
  specifying amount to fade distant points.

- shrink:

  Numeric specifying degree to which size of plotted point should
  reflect `z` position. `0` denotes no scaling; if `1`, furthest point
  will have zero size.

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
Plot3(1:10, 1:10, 1:10, cex = 7, pch = 16, axes = FALSE, asp = 1)

# Extreme values of fog and shrink will cause smallest z values to
# become invisible.
Plot3(1:10, 1:10, 1:10, cex = 7, pch = 16, axes = FALSE, asp = 1,
      fog = 1, shrink = 1)
```
