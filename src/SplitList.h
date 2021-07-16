#ifndef _TREEDIST_SPLITLIST_H
#define _TREEDIST_SPLITLIST_H

#include <stdint.h>
#include <Rcpp.h>
using namespace Rcpp;

const int16 R_BIN_SIZE = 8;

class SplitList {
public:
  int16 n_splits, n_bins;
  splitbit state[MAX_SPLITS][MAX_BINS];
  splitbit in_split[MAX_SPLITS];
  SplitList(RawMatrix);
};

#endif
