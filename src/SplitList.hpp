#ifndef _TREEDIST_SPLITLIST_H
#define _TREEDIST_SPLITLIST_H

#include <stdint.h>
#include <Rcpp.h>
#include "ints.hpp"
#include "tree_distances.hpp"

using namespace Rcpp;

const int16 R_BIN_SIZE = 8;

class SplitList {
  int16 partition(int16 lo, int16 hi);
  void swap(int16 a, int16 b);
  bool less_than(int16 a, int16 b);
  bool greater_than(int16 a, int16 b);
public:
  int16 n_splits, n_bins;
  int16 in_split[MAX_SPLITS];
  splitbit state[MAX_SPLITS][MAX_BINS];
  SplitList(RawMatrix x);
  SplitList(RawMatrix x, int16 n_tips);
  void quicksort();
  void quicksort(int16 lo, int16 hi);
};

#endif
