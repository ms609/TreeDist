# Plot a simple tree

Convenience plotting function used in vignettes and documentation.

## Usage

``` r
TreeDistPlot(
  tr,
  title = NULL,
  bold = NULL,
  leaveRoom = FALSE,
  prune = integer(0),
  graft = integer(0),
  edge.color = "black",
  edge.width = NULL,
  ...
)
```

## Arguments

- tr:

  A tree of class `phylo`, with leaves labelled as integers.

- title:

  `main` title for the plot.

- bold:

  Integer specifying which leaves to print in bold.

- leaveRoom:

  Logical specifying whether to leave space to print tree distances
  beneath the plotted tree.

- prune, graft:

  Integer vectors specifying which edges to highlight as pruned and
  grafted.

- edge.color:

  Additional parameter to `plot.phylo`; will be overridden by `prune`
  and `graft` as requested.

- edge.width, ...:

  Additional parameters to `plot.phylo`.

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)
