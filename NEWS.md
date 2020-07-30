# TreeDist 1.1.1.9000

- `TreeDistance()` and related functions now return a `dist` object when 
  computing all distances between all pairs of trees in a list.

- Set `TreeDistance()` to return a distance, rather than a similarity,
  by default (as doucmented).


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
  '[TreeSearch](https://ms609.github.io/TreeSearch)' 0.3.2.9005.
