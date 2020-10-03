# TreeDist 1.2.1.9200 (development)

- New class `ClusterTable` to allow faster distance computation with Day (1985)
  algorithm.

# TreeDist 1.2.1.9100 (development)

- Shiny app for tree space.
- `ProjectionQuality()` calculates trustworthiness and continuity of tree 
  space projections.

# TreeDist 1.2.1

- Import RdMacros package 'RdPack'.

# TreeDist 1.2.0

- `TreeDistance()` and related functions now return a `dist` object when 
  computing all distances between all pairs of trees in a list.

- Improve floating-point arithmetic in `TreeDistance()` functions.

- `TreeDistance()` now returns a distance (as documented), rather than a
  similarity.

- Fix rounding error in NNI 'Li' upper estimate, and improve NNI performance.

- Reduce precision of LAPJV so rounding errors do not result in interminable run
  times.


# TreeDist 1.1.1

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


# TreeDist 1.0.0

- Initial release, building on some draft functions included in 
  '[TreeSearch](https://ms609.github.io/TreeSearch/)' 0.3.2.9005.
