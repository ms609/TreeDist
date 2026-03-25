# Maximum Clade Information Tree

Analogous to the Maximum Clade Credibility tree: select the tree from a
posterior distribution whose clades have the highest information
content. Generate the MCC tree by specifying `info = "credibility"`.

## Usage

``` r
MCITree(trees, info = "phylogenetic", check.tips = TRUE)
```

## Arguments

- trees:

  List of `phylo` objects, optionally with class `multiPhylo`.

- info:

  Abbreviation of "phylogenetic" or "clustering", specifying the concept
  of information to employ.

- check.tips:

  Logical specifying whether to renumber leaves such that leaf numbering
  is consistent in all trees.

## Value

`MCITree()` returns the tree with the highest information content,
selected from `trees`.

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
library("TreeTools", quietly = TRUE)
trees <- as.phylo(24:40, 16)

# Maximum Clade Information tree
mci <- MCITree(trees)
SplitwiseInfo(mci)
#> [1] 133.8514
plot(mci)
p <- SplitFrequency(mci, trees) / length(trees)
LabelSplits(mci, round(p * 100), "%", bg = SupportColor(p))


# \donttest{
# Compare with Maximum Clade Credibility tree
mcc <- MCITree(trees, "credibility")
plot(mcc)
p <- SplitFrequency(mcc, trees) / length(trees)
LabelSplits(mcc, round(p * 100), "%", bg = SupportColor(p))

SplitwiseInfo(mcc)
#> [1] 128.5908
# }
```
