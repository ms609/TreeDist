# TreeDist

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![codecov](https://codecov.io/gh/ms609/TreeDist/branch/master/graph/badge.svg)](https://codecov.io/gh/ms609/TreeDist)
[![CRAN Status Badge](http://www.r-pkg.org/badges/version/TreeDist)](https://cran.r-project.org/package=TreeDist)
[![CRAN Downloads](http://cranlogs.r-pkg.org/badges/TreeDist)](https://cran.r-project.org/package=TreeDist)
[![DOI](https://zenodo.org/badge/196188301.svg)](https://zenodo.org/badge/latestdoi/196188301)

'TreeDist' is an R package that implements a suite of metrics that quantify the
topological distance between pairs of unweighted phylogenetic trees.
It also includes a simple 'Shiny' application to allow the visualization of
distance-based tree spaces, and functions to calculate the information content
of trees and splits.

'TreeDist' primarily employs metrics in the category of
'generalized Robinson–Foulds distances': they are based on comparing splits
(bipartitions) between trees, and thus reflect the relationship data within 
trees, with no reference to branch lengths.


## Generalized RF distances

The [Robinson-Foulds distance](https://ms609.github.io/TreeDist/articles/Robinson-Foulds.html)
simply tallies the number of non-trivial splits (sometimes inaccurately
termed clades, nodes or edges) that occur in both trees – any splits that are
not perfectly identical contribute one point to the distance score of zero, 
however similar or different they are.
By overlooking potential similarities between almost-identical splits, 
this conservative approach has undesirable properties.

['Generalized' RF metrics](https://ms609.github.io/TreeDist/articles/Generalized-RF.html)
generate _matchings_ that pair splits in one tree with similar splits in
the other.
Each pair of splits is assigned a similarity score; the sum of these scores in
the optimal matching then quantifies the similarity between two trees.

Different ways of calculating the the similarity between a pair of splits
lead to different tree distance metrics, implemented in the functions below:

* [`MutualClusteringInfo()`](https://ms609.github.io/TreeDist/reference/TreeDistance.html), [`SharedPhylogeneticInfo()`](https://ms609.github.io/TreeDist/reference/TreeDistance.html)
    
    Smith (2020) scores matchings based on the amount of information
    that one partition contains about the other.  The Mutual Phylogenetic
    Information assigns zero similarity to split pairs that cannot
    both exist on a single tree;  The Mutual Clustering Information metric is 
    more forgiving, and exhibits more desirable behaviour; it is the 
    recommended metric for tree comparison.
    (Its complement, 
    [`ClusteringInfoDistance()`](https://ms609.github.io/TreeDist/reference/TreeDistance.html),
    returns a tree distance.)
    
    [![Introduction to the Clustering Information Distance](man/figures/CID_talk.png)](https://durham.cloud.panopto.eu/Panopto/Pages/Viewer.aspx?id=ca5ede19-d21a-40ce-8b9e-ac6e00d7e2c0)

* [`NyeSimilarity()`](https://ms609.github.io/TreeDist/reference/NyeSimilarity.html)
    
    Nye _et al._ (2006) score matchings according to the size of the largest 
    split that is consistent with both of them, normalized against 
    the Jaccard index.  This approach is extended by B&ouml;cker _et al_. (2013)
    with the Jaccard-Robinson-Foulds metric (function 
    [`JaccardRobinsonFoulds()`](https://ms609.github.io/TreeDist/reference/JaccardRobinsonFoulds.html)).
   
* [`MatchingSplitDistance()`](https://ms609.github.io/TreeDist/reference/MatchingSplitDistance.html)
    
    Bogdanowicz and Giaro (2012) and  Lin _et al._ (2012) independently proposed
    counting the number of 'mismatched' leaves in a pair of splits.
    [`MatchingSplitInfoDistance()`](https://ms609.github.io/TreeDist/reference/TreeDistance.html)
    provides an information-based equivalent (Smith 2020).
    

The package also implements the variation of the path distance 
proposed by Kendal and Colijn (2016) (function
[`KendallColijn()`](https://ms609.github.io/TreeDist/reference/KendallColijn.html)),
approximations of the Nearest-Neighbour Interchange (NNI) distance (function
[`NNIDist()`](https://ms609.github.io/TreeDist/reference/NNIDist.html); 
following Li _et al._ (1996)), and calculates the size (function
[`MASTSize()`](https://ms609.github.io/TreeDist/reference/MASTSize.html)) and 
information content (function
[`MASTInfo()`](https://ms609.github.io/TreeDist/reference/MASTSize.html)) of the 
Maximum Agreement Subtree.

For an implementation of the Tree Bisection and Reconnection (TBR) distance, see 
the package '[TBRDist](https://ms609.github.io/TBRDist/index.html)'.

# Installation

Install and load the library from CRAN as follows:
```r
install.packages('TreeDist')
library('TreeDist')
```

You can install the development version of the package with:
```r
if (!require(devtools)) install.packages("devtools")
devtools::install_github('ms609/TreeDist')
```

# Tree space analysis

Construct tree spaces and readily visualize projected landscapes, avoiding
common analytical pitfalls (Smith, forthcoming),
using the inbuilt graphical user interface:

```r
TreeDist::MapTrees()
```

Serious analysts should consult the
[vignette](https://ms609.github.io/TreeDist/articles/treespace.html)
for a command-line interface.


# Documentation

- [Using 'TreeDist'](https://ms609.github.io/TreeDist/articles/Using-TreeDist.html)

- [Package functions](https://ms609.github.io/TreeDist/reference/index.html)

- [Tree spaces with 'TreeDist'](https://ms609.github.io/TreeDist/articles/treespace.html)

- [All vignettes](https://ms609.github.io/TreeDist/articles/)

# See also

Other R packages implementing tree distance functions include:

* '[ape](http://ape-package.ird.fr/)':
    - `cophenetic.phylo()`: Cophenetic distance
    - `dist.topo()`: Path (topological) distance, Robinson-Foulds distance.
* '[phangorn](https://cran.r-project.org/package=phangorn)'
    - `treedist()`: Path, Robinson-Foulds and approximate SPR distances.
* '[Quartet](http://ms609.github.io/Quartet/)': Triplet and Quartet distances, 
  using the tqDist algorithm.
* '[TBRDist](http://ms609.github.io/TBRDist/)': TBR and SPR distances on 
  unrooted trees, using the 'uspr' C library.
* '[distory](https://cran.r-project.org/package=distory)' (unmaintained):
  Geodesic distance

# References

- Böcker, S. _et al._ (2013) [The Generalized Robinson-Foulds
metric](https://dx.doi.org/10.1007/978-3-642-40453-5_13).
Algorithms in Bioinformatics. WABI 2013.
_Lecture Notes in Computer Science_, 8126, 156–69.

- Bogdanowicz, D. and Giaro, K. (2012) [Matching split distance for unrooted
binary phylogenetic trees](https://dx.doi.org/10.1109/TCBB.2011.48).
_IEEE/ACM Transactions on Computational Biology and Bioinformatics_, 9, 150–160. 

- Kendall, M. and Colijn, C. (2016) [Mapping phylogenetic trees to reveal
distinct patterns of evolution](https://dx.doi.org/10.1093/molbev/msw124).
_Mol Biol Evol_, 33, 2735–2743.

- Li, M., Tromp, J. and Zhang, L.-X. (1996) [Some notes on the nearest neighbour
interchange distance](https://dx.doi.org/10.1007/3-540-61332-3_168). 
_Computing and Combinatorics_, Goos, G., Hartmanis, J., Leeuwen, J., Cai, J.-Y.,
and Wong, C. K., eds. Springer, Berlin. 343–351.

- Nye, T.M.W. _et al._ (2006) [A novel algorithm and web-based tool for
comparing two alternative phylogenetic
trees](https://dx.doi.org/10.1093/bioinformatics/bti720).
_Bioinformatics_, 22, 117–119.

- Smith, M.R. (2020) [Information theoretic Generalized Robinson-Foulds
metrics for comparing phylogenetic 
trees](https://dx.doi.org/10.1093/bioinformatics/btaa614).
_Bioinformatics_, 36, 5007–5013.

- Smith, M.R. (forthcoming)
The importance of methodology when analysing landscapes of phylogenetic trees.


Please note that the 'TreeDist' project is released with a
[Contributor Code of Conduct](https://ms609.github.io/TreeDist/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
