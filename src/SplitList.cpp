#include <Rcpp.h>
using namespace Rcpp;
#include <stdint.h>
#include "tree_distances.hpp" // Before SplitList.h, for int16
#include "SplitList.hpp"

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

SplitList::SplitList(RawMatrix x, int16 n_tip) {
  n_splits = x.rows();
  const int16 
    n_input_bins = x.cols(),
    input_bins_per_bin = BIN_SIZE / R_BIN_SIZE
  ;
  
  n_bins = (n_input_bins + R_BIN_SIZE - 1) / input_bins_per_bin;
  
  if (n_bins < 1) throw std::invalid_argument("No leaves present.");
  if (n_splits < 0) throw std::invalid_argument("Negative number of splits!?");
  if (n_bins > MAX_BINS) {
    throw std::length_error("This many leaves cannot be supported. "            // # nocov
                            "Please contact the TreeDist maintainer if you "    // # nocov
                            "need to use more!");                               // # nocov
  }
  Rcout << "\n=== SplitList::SplitList(x, n_tip) ===\n";
  const int16 last_bin_tips = INLASTBIN(n_tip, R_BIN_SIZE);
  Rcout << "Last bt: " << last_bin_tips << ".\n";
  Rcout << "Last 2^bt: " << (2^last_bin_tips) << ".\n";
  const splitbit last_mask = std::pow(2, last_bin_tips) - 1;
  Rcout << "Last mask: " << last_mask << ".\n";
  
  for (int16 split = 0; split != n_splits; split++) {
    int16 last_bin = n_bins - 1;
    const int16 raggedy_bins = INLASTBIN(n_input_bins, R_BIN_SIZE);
    
    splitbit last = INSUBBIN(last_bin, 0);
    bool invert = last & splitbit(1);
    Rcout << n_input_bins << " bins in; " << raggedy_bins << " raggedy bins\n";
    Rcout << "Last bin: " << last << "; intersection with " << splitbit(1)
          << " = " << invert << "; " << (invert ? "inverting" : "Retaining")
          << ".\n";
    state[split][last_bin] = invert ? ~last & last_mask : last;
    Rcout << " State[" << split << "][" << last_bin << "] = " << state[split][last_bin] << ".\n";
    for (int16 input_bin = 1; input_bin != raggedy_bins; input_bin++) {
      Rcout << "Adding " << (splitbit (x(split, (last_bin * 2) + input_bin))) << " << "
            << (R_BIN_SIZE * input_bin) << " to state [" << split << "][" << last_bin
            << "], was " << state[split][last_bin] << "\n";
      state[split][last_bin] += invert ? 
        ~INBIN(input_bin, last_bin):
         INBIN(input_bin, last_bin);
    }
    in_split[split] = count_bits(state[split][last_bin]);
    
    for (int16 bin = 0; bin != n_bins - 1; bin++) {
      Rcout << "Split " << split << ", bin << " << bin << ".\n";
      state[split][bin] = invert ? ~INSUBBIN(bin, 0) : INSUBBIN(bin, 0);
      for (int16 input_bin = 1; input_bin != input_bins_per_bin; input_bin++) {
        Rcout << "Adding " << (splitbit (x(split, (bin * 2) + input_bin))) << " << "
              << (R_BIN_SIZE * input_bin) << " to state [" << split << "][" 
              << bin << "], was " << state[split][bin] << "\n";
        state[split][bin] += invert ?
        ~INSUBBIN(bin, input_bin):
         INSUBBIN(bin, input_bin);
      }
      in_split[split] += count_bits(state[split][bin]);
    }
  }
}
void SplitList::swap(int16 a, int16 b) {
  for (int16 bin = n_bins; bin--; ) {
    const splitbit tmp = state[a][bin];
    state[a][bin] = state[b][bin];
    state[b][bin] = tmp;
  }
}

bool SplitList::less_than(int16 a, int16 b) {
  for (int16 bin = 0; bin != n_bins; ++bin) {
    if (state[a][bin] < state[b][bin]) {
      return true;
    }
  }
  return false;
}

bool SplitList::greater_than(int16 a, int16 b) {
  for (int16 bin = 0; bin != n_bins; ++bin) {
    if (state[a][bin] > state[b][bin]) {
      return true;
    }
  }
  return false;
}

int16 SplitList::partition(int16 lo, int16 hi) {
  const int16 pivot = (hi + lo) / 2;
  for (int16 i = lo - 1, j = hi + 1; ; ) {
    do {
      ++i;
    } while (less_than(i, pivot));
    do {
      --j;
    } while (greater_than(j, pivot));
    if (i >= j) {
      return j;
    }
    swap(i, j);
  }
}

void SplitList::quicksort(int16 lo, int16 hi) {
  if (lo < hi) {
    const int16 p = partition(lo, hi);
    quicksort(lo, p);
    quicksort(p + 1, hi);
  }
}
  
void SplitList::quicksort() {
  quicksort(0, n_splits - 1);
}
