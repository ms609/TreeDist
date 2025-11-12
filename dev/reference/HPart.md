# Hierarchical partition structure

A structure of class `HPart` comprises a pointer to a C++ representation
of hierarchical partitions, with the attribute `tip.label` recording the
character labels of its leaves. `HPart` objects with identical tip
labels can be compared using
[`HierarchicalMutualInfo()`](https://ms609.github.io/TreeDist/dev/reference/HierarchicalMutualInfo.md).

## Usage

``` r
as.HPart(tree, tipLabels)

# S3 method for class 'HPart'
as.HPart(tree, tipLabels = NULL)

# Default S3 method
as.HPart(tree, tipLabels = NULL)

# S3 method for class 'list'
as.HPart(tree, tipLabels = NULL)

# S3 method for class 'phylo'
as.HPart(tree, tipLabels = TipLabels(tree))

is.HPart(x)

# S3 method for class 'HPart'
print(x, ...)

# S3 method for class 'HPart'
as.phylo(x, ...)

# S3 method for class 'HPart'
plot(x, ...)
```

## Arguments

- tree:

  An object to convert to an HPart structure, in a supported format (see
  details).

- tipLabels:

  Character vector specifying sequence in which to order tip labels.

- x:

  `HPart` object to plot.

- ...:

  Additional arguments to
  [`plot.phylo`](https://rdrr.io/pkg/ape/man/plot.phylo.html).

## Value

`HPart()` returns a structure containing a pointer to a C++
representation of a hierarchical partition structure.

## Details

An `HPart` object may be created from various representations of
hierarchical structures:

- a tree (possibly phylogenetic) of class `phylo`

- A hierarchical list of lists, in which elements are represented by
  integers 1 to n

- A vector, which will be interpreted as a flat structure in which all
  elements bearing the same label are assigned to the same cluster
