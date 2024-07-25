# TreeDist 2.8.0 (2024-07-25)

- `Islands()` allows the identification of islands of trees.

- Internal implementation of path and SPR distances, removing dependency
  on phangorn (and thus R 4.4).
  
- Add progress bar within `.MaxValue()`


# TreeDist 2.7.1 (2024-06-13)

- Documentation improvements.

- Fix `KCDiameter.multiPhylo()` for multiple trees.


# TreeDist 2.7.0 (2023-10-25)

- Fix calculation error in `StrainCol()`.

- App: Display strain in 3D tree space viewer.

- Support for distances between larger trees.

- Support unrooted trees in `VisualizeMatching()`
  ([#103](https://github.com/ms609/TreeDist/issues/103)).


# TreeDist 2.6.3 (2023-08-25)

- Fix bug when comparing a "multiPhylo" object containing a single tree.

- Documentation clarification: finding non-matching leaves.


# TreeDist 2.6.2 (2023-06-28)

- Support non-square matrices in `LAPJV()`.


# TreeDist 2.6.1 (2023-04-25)

- `StopParallel()` gains `quietly` argument to suppress unnecessary messages.

- Use "PlotTools" package for spectrum legends.

- Minor documentation tweaks.


# TreeDist 2.6.0 (2023-02-20)

- Support comparison of trees with different tips.

- Fix caching errors in `MapDist()`
  ([#98](https://github.com/ms609/TreeDist/issues/98)).

- Update tests for compatibility with ape 5.7.


# TreeDist 2.5.0 (2022-10-07)

- New functions to measure cluster sizes (see [`?"cluster-statistics"`](
  https://ms609.github.io/TreeDist/reference/cluster-statistics.html)).

- `KMeansPP()` conducts clustering using K-means++, replacing K-means in app.

- New [vignette](https://ms609.github.io/TreeDist/articles/landscapes.html)
  on tree landscape analysis.


# TreeDist 2.4.1 (2022-07-20)

- New [vignette](https://ms609.github.io/TreeDist/articles/compare-treesets.html)
  on how to compare tree sets.

- `PathVector()` now treats trees with a root node as rooted.

- Fix plot layout in [treespace vignette](https://ms609.github.io/TreeDist/articles/treespace.html).

- Informative failure when not enough memory for `consensus_info()`.

- Replace `throw` with `stop` in C++.


# TreeDist 2.4.0 (2022-03-23)

- Correct calculation of trustworthiness and continuity metrics.

- Depict strain in minimum spanning trees with `StrainCol()` and helper
  function `MSTSegments()`.

- Update tests for consistency with "TreeTools" v1.7.

- Use lighter Rcpp headers.


# TreeDist 2.3.0 (2022-01-04)

- Support `ConsensusInfo(p > 0.5)`.

- Address hypervolume comparison in vignettes.

- Support uniform manifold approximation and projection in app.


# TreeDist 2.2.0 (2021-09-13)

- Speed improvements, using optimizations suggested by Alexis Stamatakis'
  Bioinformatics group.
  
- Support for parallel computation via `StartParallel()`.
  
- Progress bars.


# TreeDist 2.1.1 (2021-07-13)

- Solaris compatibility.

- Modest vignette improvements.

- spic/scic abbreviation recognition.


# TreeDist 2.1.0 (2021-07-12)

## New features

- `ConsensusInfo()` quickly calculates the splitwise information content of the
  consensus of a set of trees, after Smith (forthcoming).

- `SplitwiseInfo()` and `ClusteringInfo()` gain a `p` parameter to reflect the
  reduced information content of splits with lower support values, and a `sum`
  parameter to allow return of individual split information content.

- `KCDiameter()` approximates the diameter of the Kendall-Colijn metric.

- `Plot3()` (experimental) provides pseudo-3D plotting.


## Renamed functions

- `Project()`/`ProjectionQuality()` re-named to `MapTrees()`/`MappingQuality()`.

- `SpectralClustering()` re-named to `SpectralEigens()`.


## Improvements

- Add self-organizing map example to treespace vignette.

- Allow the specification of custom vectors in the Kendall--Colijn metric.

- Faster all-to-all tree distance calculation.


# TreeDist 2.0.3 (2021-01-31)

- Diagnose and fix memory leaks, including over-long reported matchings.

- Explicitly import shiny/shinyjs functions.


# TreeDist 2.0.0 (2021-01-20)

- `Project()` launches 'shiny' app for projection and analysis of tree space.

- `ProjectionQuality()` calculates trustworthiness and continuity of tree 
  space mappings.
  
- Faster calculation of Robinson–Foulds distance (using algorithm of Day (1985))
  and clustering information distance.
  
- New class `ClusterTable` to allow faster distance computation with Day (1985)
  algorithm.
  
- Improve error messages in `CalculateTreeDist()`.

- Improvements to vignettes.

- Use package 'vdiffr' conditionally.


# TreeDist 1.2.1 (2020-09-17)

- Import RdMacros package 'RdPack'.


# TreeDist 1.2.0 (2020-08-28)

- `TreeDistance()` and related functions now return a `dist` object when 
  computing all distances between all pairs of trees in a list.

- Improve floating-point arithmetic in `TreeDistance()` functions.

- `TreeDistance()` now returns a distance (as documented), rather than a
  similarity.

- Fix rounding error in NNI 'Li' upper estimate, and improve NNI performance.

- Reduce precision of LAPJV so rounding errors do not result in interminable run
  times.


# TreeDist 1.1.1 (2020-07-10)

- Fix range errors when calculating tree distances.


# TreeDist 1.1.0

- Improvements to `NNIDist()` in light of Fack _et al._ (2002).

- Add `NNIDiameter()`: approximate diameter of NNI distance.
 
- Remove vignette 'Interpreting tree distances': duplicates
  https://ms609.github.io/TreeDistData/articles/09-expected-similarity.html.
  
- Remove redundant data object `oneOverlap`.

- Fix an issue when installing on R 3.x (require C++11 to ensure declaration 
  of `UINT_FAST16_MAX`).
  
- Fix memory-handling bug in `lapjv()`.


# TreeDist 1.0.0 (2020-06-30)

- Initial release, building on some draft functions included in 
  '[TreeSearch](https://ms609.github.io/TreeSearch/)' 0.3.2.9005.
