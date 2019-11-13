#include <Rcpp.h>
using namespace Rcpp;
#include <stdint.h>

const int R_BIN_SIZE = 8;

class SplitList {
public:
  int n_splits, n_bins;
  splitbit state[MAX_SPLITS][MAX_BINS];
  SplitList(NumericMatrix);
};
