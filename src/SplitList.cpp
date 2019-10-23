#include <Rcpp.h>
using namespace Rcpp;
#include <stdint.h>
#include "SplitList.h"

__attribute__((constructor))
  void initialize_array()
  {
    for (int32_t i = 0; i < 65536; i++) {
      int32_t n_bits = 0;
      for (int j = 0; j < 16; j++) {
        if ((i & powers_of_two[j])) ++n_bits;
      }
      bitcounts[i] = n_bits;
    }
  }

int count_bits_32 (uint32_t x) {
  return bitcounts[x & right16bits] + bitcounts[x >> 16];
}

SplitList::SplitList(NumericMatrix x) {
  n_splits = x.rows();
  n_bins = x.cols();
  
  if (n_bins < 1) throw std::invalid_argument("No tips present.");
  if (n_splits < 1) throw std::invalid_argument("No splits present.");
  if (n_bins > 100) {
    throw std::length_error("No more than 3200 tips can be supported. Please contact the TreeDist maintainer if you need to use more!");
  }
  
  for (int i = 0; i < n_splits; i++) {
    for (int j = 0; j < n_bins; j++) {
      state[i][j] = (uint32_t) x(i, j);
    }
  }
}