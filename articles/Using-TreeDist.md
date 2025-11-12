# Calculate tree similarity with 'TreeDist'

This document should contain all you need to get started measuring tree
distances with ‘TreeDist’. If you get stuck, please [let me
know](https://github.com/ms609/TreeDist/issues/new?title=Suggestion:+)
so I can improve this documentation.

## Loading trees

Instructions for loading phylogenetic trees into R can be found in a
[separate
vignette](https://ms609.github.io/TreeTools/articles/load-trees.html).
For these examples, we’ll enter two simple trees by hand:

``` r
tree1 <- ape::read.tree(text = '(A, ((B, (C, (D, E))), ((F, G), (H, I))));')
tree2 <- ape::read.tree(text = '(A, ((B, (C, (D, (H, I)))), ((F, G), E)));')
```

## Calculating distances

We can calculate distances between pairs of trees using the ‘TreeDist’
package.

First we’ll install the package. We can either install the stable
version from the CRAN repository:

``` r
install.packages('TreeDist')
```

or the development version, from GitHub – which will contain the latest
features but may not be as extensively tested:

``` r
devtools::install_github('ms609/TreeDist')
```

Then we’ll load the package in to R’s working environment:

``` r
library('TreeDist')
```

Now the package’s functions are available within R. Let’s proceed to
calculate some tree distances.

### Pairs of trees

Calculating the distance between two trees is as simple as:

``` r
distance <- TreeDistance(tree1, tree2)
```

The convenience function
[`TreeDistance()`](https://ms609.github.io/TreeDist/reference/TreeDistance.md)
returns the variation of clustering information between two trees,
[normalized](https://ms609.github.io/TreeDist/articles/using-distances.html#normalizing)
against the total information content of all splits.

### Multiple comparisons

If you have more than two trees to compare, you can send a list of trees
(class: `list` or `multiPhylo`) to the distance comparison function. The
function will then calculate the distance between each tree in the first
list and each tree in the second.

``` r
oneTree <- ape::rtree(11)
twoTrees <- structure(list(one = ape::rtree(11), two = ape::rtree(11)),
                      class = 'multiPhylo')
threeTrees <- list(a = ape::rtree(11), b = ape::rtree(11), c = ape::rtree(11))

TreeDistance(oneTree, twoTrees)
```

    ##       one       two 
    ## 0.7826864 0.6856279

``` r
TreeDistance(twoTrees, threeTrees)
```

    ##             a         b         c
    ## one 0.8620894 0.7514794 0.7861809
    ## two 0.8556624 0.7150163 0.7658788

## Visualizing a matching

[Generalized Robinson–Foulds
metrics](https://ms609.github.io/TreeDist/articles/Generalized-RF.md),
such as the variation of clustering information, rely on matching each
split within a tree with another split in the other tree.  
We can view an optimal matching:

``` r
VisualizeMatching(ClusteringInfoDistance, tree1, tree2)
```

![Pair of similar phylogenetic trees with matched splits highlighted
according to the amount of clustering information in
common.](Using-TreeDist_files/figure-html/visualise-matching-1.png)

This shows the six splits in tree 1, and the paired splits in tree
two.  
Each split is labelled with a measure of its similarity, which is its
contribution to the total tree similarity score.

We can view this information in a format accessible for further
examination in R with:

``` r
ClusteringInfoDistance(tree1, tree2, reportMatching = TRUE)
```

    ## [1] 6.960578
    ## attr(,"matching")
    ## [1] 1 2 3 5 6 4
    ## attr(,"matchedSplits")
    ## [1] "B C D E | A F G H I => B C D H I | A E F G"
    ## [2] "C D E | A B F G H I => C D H I | A B E F G"
    ## [3] "D E | A B C F G H I => D H I | A B C E F G"
    ## [4] "F G H I | A B C D E => E F G | A B C D H I"
    ## [5] "F G | A B C D E H I => F G | A B C D E H I"
    ## [6] "H I | A B C D E F G => H I | A B C D E F G"
    ## attr(,"matchedScores")
    ## [1] 0.09109101 0.07278023 0.02475761 0.07278023 0.76420451 0.76420451
    ## attr(,"pairScores")
    ##             [,1]        [,2]       [,3]       [,4]       [,5]       [,6]
    ## [1,] 0.091091008 0.007214618 0.01831078 0.22478751 0.01831078 0.22478751
    ## [2,] 0.018310782 0.072780226 0.00000000 0.15200728 0.00000000 0.15200728
    ## [3,] 0.002565287 0.002565287 0.02475761 0.09288851 0.02475761 0.09288851
    ## [4,] 0.007214618 0.007214618 0.07278023 0.31976006 0.07278023 0.31976006
    ## [5,] 0.319760062 0.224787510 0.15200728 0.09288851 0.45810590 0.76420451
    ## [6,] 0.224787510 0.319760062 0.45810590 0.76420451 0.15200728 0.09288851

Here, the `pairScores` attribute lists the score of each possible
matching of splits.

We can identify the splits with:

``` r
splits <- as.character(TreeTools::as.Splits(tree2))
splits
```

    ##                    12                    13                    14 
    ## "B C D H I | A F G E" "C D H I | A B F G E" "D H I | A B C F G E" 
    ##                    15                    16                    17 
    ## "H I | A B C D F G E" "F G E | A B C D H I" "F G | A B C D H I E"

The names of the splits correspond to the number of an associated node
in the original tree:

``` r
oldPar <- par(mar = rep(0, 4))
plot(tree2)
ape::nodelabels()
ape::nodelabels(splits, as.integer(names(splits)), 
                adj = c(1.1, -0.2), cex = 0.8, frame = 'none')
```

![Phylogenetic tree with nodes numbered, and labelled with the splits to
which they
correspond.](Using-TreeDist_files/figure-html/named-splits-1.png)

Note that strictly, (informative) splits are associated with (internal)
edges. To avoid listing the same split twice, nodes close to the root
(here, 10 and 11) will not be associated with a split.

## What next?

You may wish to:

- [Provide
  context](https://ms609.github.io/TreeDist/articles/using-distances.md)
  for tree distances

- Compare trees with [different
  tips](https://ms609.github.io/TreeDist/articles/different-leaves.md)

- Review [available distance
  measures](https://ms609.github.io/TreeDist/index.html) and the
  corresponding
  [functions](https://ms609.github.io/TreeDist/reference/index.html#section-tree-distance-measures)

- [Interpret tree distance
  metrics](https://ms609.github.io/TreeDistData/articles/09-expected-similarity.html)

- Visualize [tree
  landscapes](https://ms609.github.io/TreeDist/articles/treespace.md)
  using distance-based tree spaces
