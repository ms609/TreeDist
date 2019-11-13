#include <Rcpp.h>
using namespace Rcpp;
#include "tree_distances.h"
#include <stdint.h>
#include "SplitList.h"

SplitList::SplitList(NumericMatrix x) {
  n_splits = x.rows();
  const int n_input_bins = x.cols(),
    bin_ratio = BIN_SIZE / R_BIN_SIZE;
  n_bins = (n_input_bins + R_BIN_SIZE - 1) / bin_ratio;
  
  if (n_bins < 1) throw std::invalid_argument("No tips present.");
  if (n_splits < 1) throw std::invalid_argument("No splits present.");
  if (n_bins > MAX_BINS * 2) {
    throw std::length_error("This many tips cannot be supported. Please contact the TreeDist maintainer if you need to use more!");
    /*throw std::length_error(printf("No more than %i tips can be supported. Please contact the TreeDist maintainer if you need to use more!",
                                   MAX_TIPS));*/
  }
  
  for (int i = 0; i != n_splits; i++) {
    for (int j = 0; j != n_bins - 1; j++) {
      state[i][j] = (splitbit) x(i, (j * 2));
      for (int input_bin = 1; input_bin != bin_ratio; input_bin++) {
        state[i][j] += ((splitbit (x(i, (j * 2) + input_bin)))
                           << (R_BIN_SIZE * input_bin));
      }
    }
    int j = n_bins - 1;
    state[i][j] = x(i, j * 2);
    const int raggedy_bins = R_BIN_SIZE - ((R_BIN_SIZE - n_input_bins) % R_BIN_SIZE);
    for (int input_bin = 1; input_bin != raggedy_bins; input_bin++) {
      state[i][j] += ((splitbit (x(i, (j * 2) + input_bin)))
                        << (R_BIN_SIZE * input_bin));
    }
  }
}