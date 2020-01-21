---
title: "Calculate tree similarity with TreeDist"
author: "Martin R. Smith"
date: "2020-01-21"
output: rmarkdown::html_vignette
bibliography: ../inst/REFERENCES.bib
csl: ../inst/apa-old-doi-prefix.csl
vignette: >
  %\VignetteIndexEntry{Using TreeDist}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

This document should contain all you need to get started measuring tree 
distances with TreeDist.  If you get stuck, please 
[let me know](https://github.com/ms609/TreeDist/issues/new?title=Suggestion:+)
so I can improve this documentation.

## Loading trees

Instructions for loading phylogenetic trees into R can be found in a [separate vignette](https://ms609.github.io/TreeTools/articles/load-trees.html).
For these examples, we'll enter two simple trees by hand:


```r
tree1 <- ape::read.tree(text='(A, ((B, (C, (D, E))), ((F, G), (H, I))));')
tree2 <- ape::read.tree(text='(A, ((B, (C, (D, (H, I)))), ((F, G), E)));')
```

## Calculating distances

We can calculate distances between pairs of trees using the `TreeDist` package.

First we'll install the package.
We can either install the stable version from the CRAN repository:

```r
install.packages('TreeDist')
```

or the development version, from GitHub -- which will contain the latest
features but may not be as extensively tested:

```r 
devtools::install_github('ms609/TreeDist')
```

Then we'll load the package in to R's working environment:


```r
library('TreeDist')
```

Now the package's functions are available within R.
Let's proceed to calculate some tree distances.

### Pairs of trees

Calculating the distance between two trees is as simple as:


```r
distance <- TreeDistance(tree1, tree2)
```

The convenience function `TreeDistance` returns the variation of clustering
information between two trees, [normalized](using-distances.html)
against the total information content of all splits.  
For more on different tree distance metrics, see the
[tree distance vignette](tree-distances.html).

### Multiple comparisons

If you have more than two trees to compare, you can send a list of trees 
(class: `list` or `multiPhylo`) to the distance comparison function.
The function will then calculate the distance between each tree in the first
list and each tree in the second.


```r
oneTree <- ape::rtree(11)
twoTrees <- structure(list(one = ape::rtree(11), two = ape::rtree(11)),
                      class='multiPhylo')
threeTrees <- list(a = ape::rtree(11), b = ape::rtree(11), c = ape::rtree(11))

TreeDistance(oneTree, twoTrees)
```

```
##       one       two 
## 0.1314848 0.5287732
```

```r
TreeDistance(twoTrees, threeTrees)
```

```
##             a         b         c
## one 0.2256777 0.2277760 0.1641394
## two 0.1989870 0.3299788 0.1819315
```

## Visualizing a matching

[Generalized Robinson-Foulds metrics](Generalized-RF.html),
such as the variation of clustering information,
rely on matching each split within a tree with another split in the other tree.  
We can view an optimal matching:


```r
VisualizeMatching(ClusteringInfoDistance, tree1, tree2)
```

<img src="C:/Users/ms609/AppData/Local/Temp/Rtmpia1ydX/preview-c2c55264837.dir/Using-TreeDist_files/figure-html/visualise-matching-1.png" width="90%" style="display: block; margin: auto;" />

This shows the eight splits in tree 1, and the paired splits in tree two.  
Each split is labelled with a measure of its similarity, which is its contribution 
to the total tree similarity score.

We can view this information in a format accessible for further examination in R with:


```r
ClusteringInfoDistance(tree1, tree2, reportMatching = TRUE)
```

```
## [1] 62.6452
## attr(,"matching")
## [1] 1 2 3 5 6 4
## attr(,"pairScores")
##            [,1]       [,2]      [,3]      [,4]      [,5]      [,6]
## [1,] 0.81981907 0.06493157 0.1647970 2.0230876 0.1647970 2.0230876
## [2,] 0.16479704 0.65502203 0.0000000 1.3680656 0.0000000 1.3680656
## [3,] 0.02308759 0.02308759 0.2228185 0.8359966 0.2228185 0.8359966
## [4,] 0.06493157 0.06493157 0.6550220 2.8778406 0.6550220 2.8778406
## [5,] 2.87784056 2.02308759 1.3680656 0.8359966 4.1229531 6.8778406
## [6,] 2.02308759 2.87784056 4.1229531 6.8778406 1.3680656 0.8359966
## attr(,"matchedSplits")
## [1] "B C D E | A F G H I => B C D H I | A E F G"
## [2] "C D E | A B F G H I => C D H I | A B E F G"
## [3] "D E | A B C F G H I => D H I | A B C E F G"
## [4] "F G H I | A B C D E => E F G | A B C D H I"
## [5] "F G | A B C D E H I => F G | A B C D E H I"
## [6] "H I | A B C D E F G => H I | A B C D E F G"
```

Here, the `pairScores` attribute lists the score of each possible matching of splits.

We can identify the splits with:


```r
splits <- as.character(TreeTools::as.Splits(tree2))
splits
```

```
##                    12                    13                    14 
## "B C D H I | A F G E" "C D H I | A B F G E" "D H I | A B C F G E" 
##                    15                    16                    17 
## "H I | A B C D F G E" "F G E | A B C D H I" "F G | A B C D H I E"
```

The names of the splits correspond to the number of the nodes in the original tree:


```r
par(mar = rep(0, 4))
plot(tree2)
ape::nodelabels()
ape::nodelabels(splits, as.integer(names(splits)), 
                adj=c(1.1, -0.2), cex=0.8, frame='none')
```

<img src="C:/Users/ms609/AppData/Local/Temp/Rtmpia1ydX/preview-c2c55264837.dir/Using-TreeDist_files/figure-html/named-splits-1.png" width="80%" style="display: block; margin: auto;" />

## What next?

You may wish to:

- [Provide context](using-distances.html) for tree distances

- Understand [how different metrics work](tree-distances.html)
