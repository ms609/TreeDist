#include <Rcpp.h>
using namespace Rcpp;
#include "tree_distances.h"
#include <stdint.h>
#include "SplitList.h"

SplitList::SplitList(RawMatrix x) {
  n_splits = x.rows();
  const int n_input_bins = x.cols(),
    input_bins_per_bin = BIN_SIZE / R_BIN_SIZE;
  
  n_bins = (n_input_bins + R_BIN_SIZE - 1) / input_bins_per_bin;
  
  if (n_bins < 1) throw std::invalid_argument("No tips present.");
  if (n_splits < 0) throw std::invalid_argument("Negative number of splits!?");
  if (n_bins > MAX_BINS) {
    throw std::length_error("This many tips cannot be supported. "
                            "Please contact the TreeDist maintainer if you "
                            "need to use more!");
  }
  
  for (int split = 0; split != n_splits; split++) {
    for (int bin = 0; bin != n_bins - 1; bin++) {
      /*Rcout << "Split " << split << ", bin << " << bin << ".\n";*/
      state[split][bin] = (splitbit) x(split, (bin * input_bins_per_bin));
      for (int input_bin = 1; input_bin != input_bins_per_bin; input_bin++) {
        /*Rcout << "Adding " << (splitbit (x(split, (bin * 2) + input_bin))) << " << "
              << (R_BIN_SIZE * input_bin) << " to state [" << split << "][" 
              << bin << "], was " << state[split][bin] << "\n";*/
        state[split][bin] += ((splitbit (x(split, (bin * input_bins_per_bin) + input_bin)))
                                << (R_BIN_SIZE * input_bin));
      }
    }
    int bin = n_bins - 1;
    const int raggedy_bins = R_BIN_SIZE - 
      ((R_BIN_SIZE - (n_input_bins % R_BIN_SIZE)) % R_BIN_SIZE);
    /*Rcout << n_input_bins << " bins in; " << raggedy_bins << " raggedy bins\n";*/
    state[split][bin] = x(split, bin * input_bins_per_bin);
    /*Rcout << " State[" << split << "][" << bin << "] = " << state[split][bin] << ".\n";*/
    for (int input_bin = 1; input_bin != raggedy_bins; input_bin++) {
      /*Rcout << "Adding " << (splitbit (x(split, (bin * 2) + input_bin))) << " << "
            << (R_BIN_SIZE * input_bin) << " to state [" << split << "][" << bin 
            << "], was " << state[split][bin] << "\n";*/
      state[split][bin] += ((splitbit (x(split, (bin * input_bins_per_bin) + input_bin)))
                        << (R_BIN_SIZE * input_bin));
    }
  }
}

SplitList::SplitList(std::vector<int> x, const int n_input_bins) {
  n_splits = x.size() / n_input_bins;
  const int input_bins_per_bin = BIN_SIZE / R_BIN_SIZE;
  n_bins = (n_input_bins + R_BIN_SIZE - 1) / input_bins_per_bin;
  
  if (n_bins < 1) throw std::invalid_argument("No tips present.");
  if (n_splits < 1) throw std::invalid_argument("No splits present.");
  if (n_bins > MAX_BINS) {
    throw std::length_error("This many tips cannot be supported. Please contact the TreeDist maintainer if you need to use more!");
  }
  
  for (int split = 0; split != n_splits; split++) {
    for (int bin = 0; bin != n_bins - 1; bin++) {
      state[split][bin] = (splitbit) x[split + (n_splits * bin * input_bins_per_bin)];
      for (int input_bin = 1; input_bin != input_bins_per_bin; input_bin++) {
        state[split][bin] += ((splitbit (x[(split + (n_splits * (bin * input_bins_per_bin)) + input_bin)]))
                                << (R_BIN_SIZE * input_bin));
      }
    }
    int bin = n_bins - 1;
    const int raggedy_bins = R_BIN_SIZE - ((R_BIN_SIZE - (n_input_bins % R_BIN_SIZE)) % R_BIN_SIZE);
    state[split][bin] = x[split + (n_splits * bin * input_bins_per_bin)];
    for (int input_bin = 1; input_bin != raggedy_bins; input_bin++) {
      state[split][bin] += ((splitbit (x[split + (n_splits * ((bin * input_bins_per_bin)) + input_bin)]))
                        << (R_BIN_SIZE * input_bin));
    }
  }
}