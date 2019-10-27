#include <Rcpp.h>
using namespace Rcpp;
#include "tree_distances.h"
#include <stdint.h>
#include "SplitList.h"

SplitList::SplitList(NumericMatrix x) {
  n_splits = x.rows();
  n_bins = x.cols();
  
  if (n_bins < 1) throw std::invalid_argument("No tips present.");
  if (n_splits < 1) throw std::invalid_argument("No splits present.");
  if (n_bins > MAX_BINS) {
    throw std::length_error("This many tips cannot be supported. Please contact the TreeDist maintainer if you need to use more!");
    /*throw std::length_error(printf("No more than %i tips can be supported. Please contact the TreeDist maintainer if you need to use more!",
                                   MAX_TIPS));*/
  }
  
  for (int i = 0; i < n_splits; i++) {
    for (int j = 0; j < n_bins; j++) {
      state[i][j] = (splitbit) x(i, j);
    }
  }
}