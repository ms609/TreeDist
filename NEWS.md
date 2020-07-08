# TreeDist 1.1.0.9000 (development)

- Fix integer range error when calculating tree distances.

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
