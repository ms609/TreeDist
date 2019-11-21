#include <Rcpp.h>
#include "tree_distances.h"
#include "splits.h"
#include "SplitList.h"

using namespace Rcpp;

// [[Rcpp::export]]
List cpp_nni_distance (IntegerMatrix edge1, IntegerMatrix edge2,
                       IntegerVector nTip) {
  RawMatrix splits1 = cpp_edge_to_splits(edge1, nTip);
  RawMatrix splits2 = cpp_edge_to_splits(edge2, nTip);
  
  IntegerVector match = cpp_robinson_foulds_distance(splits1, splits2, nTip)(2);
  
  
  List ret = List::create(match);
  return (ret);
}