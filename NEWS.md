# TreeDist 1.0.0.9001 (development)

 - Improvements to `NNIDist()` in light of Fack _et al._ (2002).

 - Add `NNIDiameter()`: approximate diameter of NNI distance.

 - Fix an issue when installing on R 3.x (require C++11 to ensure declaration 
   of `UINT_FAST16_MAX`).
   
 - Fix memory-handling bug in `lapjv()`.
  
 - Remove vignette 'Interpreting tree distances': duplicates
   https://ms609.github.io/TreeDistData/articles/09-expected-similarity.html.
- Remove redundant data object `oneOverlap`.

# TreeDist 1.0.0

 - Initial release, building on some draft functions included in 
   '[TreeSearch](https://ms609.github.io/TreeSearch)' 0.3.2.9005.