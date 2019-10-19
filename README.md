# TreeDist
<!--#TODO: https://devguide.ropensci.org/building.html#readme -->

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)<!--[![Project Status: Inactive – The project has reached a stable, usable state but is no longer being actively developed; support/maintenance will be provided as time allows.](http://www.repostatus.org/badges/latest/inactive.svg)](http://www.repostatus.org/#inactive)-->
[![Travis build status](https://travis-ci.org/ms609/TreeDist.svg?branch=master)](https://travis-ci.org/ms609/TreeDist)
[![codecov](https://codecov.io/gh/ms609/TreeDist/branch/master/graph/badge.svg)](https://codecov.io/gh/ms609/TreeDist)
[![CRAN Status Badge](http://www.r-pkg.org/badges/version/TreeDist)](https://cran.r-project.org/package=TreeDist)
[![CRAN Downloads](http://cranlogs.r-pkg.org/badges/TreeDist)](https://cran.r-project.org/package=TreeDist)
<!--[![DOI](https://zenodo.org/badge/98171642.svg)](https://zenodo.org/badge/latestdoi/98171642)-->

`TreeDist` implements a suite of metrics that quantify the topological 
distance between pairs of unweighted phylogenetic trees.
The metrics generally fall in the category of "generalized Robinson-Foulds
distances": they are based on comparing partitions between trees, and thus
reflect the relationship data within trees, with no reference to branch lengths.

The [Robinson-Foulds distance](https://ms609.github.io/TreeDist/articles/Robinson-Foulds.html) simply tallies the number of bipartition splits
(loosely, clades) that occur in both trees -- any clades that are not perfectly 
identical are assigned a score of zero, however similar or different they are.
By overlooking potential similarities between almost-identical bipartitions, 
this conservative approach has undesirable properties.

['Generalized' RF metrics](https://ms609.github.io/TreeDist/articles/Generalized-RF.html)
pair each bipartition in one tree with a similar
bipartition in the other.  Each pair of partitions is assigned a similarity 
score; the sum of these scores in the optimal partition matching 
then describes the similarity between two trees. 

Different ways of calculating the the similarity between a pair of partitions 
lead to different tree distance metrics, implemented in the functions below:


* [`MutualClusteringInfo`](https://ms609.github.io/TreeDist/reference/TreeDistance.html), [`MutualPhylogeneticInfo`](https://ms609.github.io/TreeDist/reference/TreeDistance.html)
    
    Smith (forthcoming) scores matchings based on the amount of information
    that one partition contains about the other.  The Mutual Phylogenetic
    Information imposes arboreal matching: i.e. a pair of splits that cannot
    both exist on a single tree are assigned zero similarity.  The Mutual 
    Clustering Information metric is more forgiving, and exhibits more 
    desirable behaviour; it is the recommended metric for tree comparison.
    (Its complement, [`VariationOfClusteringInfo`](https://ms609.github.io/TreeDist/reference/TreeDistance.html), returns a tree 
    distance.)

* [`NyeTreeSimilarity`](https://ms609.github.io/TreeDist/reference/NyeTreeSimilarity.html)
    
    Nye _et al._ (2006) score matchings according to the size of the largest 
    bipartition that is consistent with both of them, normalized against 
    the Jaccard index.  This approach is extended by B&ouml;cker _et al_. (2013)
    with the Jaccard-Robinson-Foulds metric (function 
    [`JaccardRobinsonFoulds`](https://ms609.github.io/TreeDist/reference/JaccardRobinsonFoulds.html)).
   
* [`MatchingSplitDistance`](https://ms609.github.io/TreeDist/reference/MatchingSplitDistance.html)
    
    Bogdanowicz and Giaro (2012), and independently Lin _et al._ (2012), count 
    the number of 'mismatched' terminals in a pair of bipartitions.

The package also implements the variation of the path distance 
proposed by Kendal and Colijn (2016) (function [`KendallColijn`](https://ms609.github.io/TreeDist/reference/KendallColijn.html)).

# Installation

<!--
#TODO: submit to CRAN!
Install and load the library from CRAN as follows:
```r
install.packages('TreeDist')
library('TreeDist')
```

If you're feeling brave, y-->You can install the development version of the package with:
```r
if (!require(devtools)) install.packages("devtools")
devtools::install_github('ms609/TreeDist')
```
# Documentation

- [Using TreeDist](https://ms609.github.io/TreeDist/articles/Using-TreeDist.html)

- [Package functions](https://ms609.github.io/TreeDist/reference/index.html)

- [All vignettes](https://ms609.github.io/TreeDist/articles/)

# See also

Other tree distance functions are implemented in:

* [ape](http://ape-package.ird.fr/):
    - `cophenetic.phylo`: Cophenetic distance
    - `dist.topo`: Path (topological) distance, Robinson-Foulds distance.
* [phangorn](https://cran.r-project.org/package=phangorn)
    - `treedist`: SPR, Robinson-Foulds and path distances.
* [Quartet](http://ms609.github.io/Quartet/): Quartet and Triplet distances, 
  using the tqDist algorithm.
* [distory](https://cran.r-project.org/package=distory) (unmaintained): Geodistic distance

# References

- Böcker, S. _et al._ (2013) [The Generalized Robinson-Foulds
Metric](https://dx.doi.org/10.1007/978-3-642-40453-5_13). Algorithms in Bioinformatics. WABI 2013. Lecture Notes in Computer Science, 8126, 156–69.

- Bogdanowicz, D. and Giaro, K. (2012) [Matching split distance for unrooted
binary phylogenetic trees](https://dx.doi.org/10.1109/TCBB.2011.48). IEEE/ACM Transactions on Computational Biology and Bioinformatics, 9, 150–160. 

- Kendall, M. and Colijn, C. (2016) [Mapping phylogenetic trees to reveal
distinct patterns of evolution](https://dx.doi.org/10.1093/molbev/msw124). Mol Biol Evol, 33, 2735–2743.

- Nye, T.M.W. _et al._ (2006) [A novel algorithm and web-based tool for
comparing two alternative phylogenetic trees](https://dx.doi.org/10.1093/bioinformatics/bti720). Bioinformatics, 22, 117–119.

- Smith, M.R. (in review) [Information theoretic Generalized Robinson-Foulds
metrics for comparing phylogenetic trees].

Please note that the 'TreeDist' project is released with a
[Contributor Code of Conduct](CODE_OF_CONDUCT.md).
By contributing to this project, you agree to abide by its terms.
