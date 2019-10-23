#include <Rcpp.h>
using namespace Rcpp;
#include <stdint.h>

const uint32_t right16bits = 65535;
const uint32_t powers_of_two[16] = {1, 2, 4, 8, 16, 32, 64, 128, 256, 512,
                                    1024, 2048, 4096, 8192, 16384, 32768};

int count_bits_32 (uint32_t);

class SplitList {
public:
  int n_splits, n_bins;
  uint32_t state[3200][100]; /* Maximum tips supported: 3200 */
  SplitList(NumericMatrix);
};
