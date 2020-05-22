#include <stdint.h>
#include <Rcpp.h>
using namespace Rcpp;

const int16 R_BIN_SIZE = 8;

class SplitList {
public:
  int16 n_splits, n_bins;
  splitbit state[MAX_SPLITS][MAX_BINS];
  SplitList(RawMatrix);
  SplitList(std::vector<int16>, const int16);
};
