#include <Rcpp.h>
using namespace Rcpp;
#include "tree_distances.h" // Before SplitList.h, for int16
#include <stdint.h>
#include "SplitList.h"

SplitList::SplitList(RawMatrix x) {
  n_splits = x.rows();
  const int16 n_input_bins = x.cols(),
              input_bins_per_bin = BIN_SIZE / R_BIN_SIZE;

  n_bins = (n_input_bins + R_BIN_SIZE - 1) / input_bins_per_bin;
  
  if (n_bins < 1) throw std::invalid_argument("No tips present.");
  if (n_splits < 0) throw std::invalid_argument("Negative number of splits!?");
  if (n_bins > MAX_BINS) {
    throw std::length_error("This many tips cannot be supported. "              // # nocov
                            "Please contact the TreeDist maintainer if you "    // # nocov
                            "need to use more!");                               // # nocov
  }
  
  for (int16 split = 0; split != n_splits; split++) {
    int16 last_bin = n_bins - 1;
    const int16 raggedy_bins = R_BIN_SIZE - 
      ((R_BIN_SIZE - (n_input_bins % R_BIN_SIZE)) % R_BIN_SIZE);
    /*Rcout << n_input_bins << " bins in; " << raggedy_bins << " raggedy bins\n";*/
    state[split][last_bin] = x(split, last_bin * input_bins_per_bin);
    /*Rcout << " State[" << split << "][" << bin << "] = " << state[split][bin] << ".\n";*/
    for (int16 input_bin = 1; input_bin != raggedy_bins; input_bin++) {
      /*Rcout << "Adding " << (splitbit (x(split, (bin * 2) + input_bin))) << " << "
              << (R_BIN_SIZE * input_bin) << " to state [" << split << "][" << bin 
              << "], was " << state[split][bin] << "\n";*/
      state[split][last_bin] += ((splitbit (x(split, (last_bin * input_bins_per_bin) + input_bin)))
                                   << (R_BIN_SIZE * input_bin));
    }
    in_split[split] = count_bits(state[split][last_bin]);
    
    for (int16 bin = 0; bin != n_bins - 1; bin++) {
      /*Rcout << "Split " << split << ", bin << " << bin << ".\n";*/
      state[split][bin] = (splitbit) x(split, (bin * input_bins_per_bin));
      for (int16 input_bin = 1; input_bin != input_bins_per_bin; input_bin++) {
        /*Rcout << "Adding " << (splitbit (x(split, (bin * 2) + input_bin))) << " << "
              << (R_BIN_SIZE * input_bin) << " to state [" << split << "][" 
              << bin << "], was " << state[split][bin] << "\n";*/
        state[split][bin] += ((splitbit (x(split, (bin * input_bins_per_bin) + input_bin)))
                                << (R_BIN_SIZE * input_bin));
      }
      in_split[split] += count_bits(state[split][bin]);
    }
  }
}
