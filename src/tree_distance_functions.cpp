#include <Rcpp/Lightest>
#include <TreeTools/SplitList.h> /* for SL_MAX_TIPS */

#include <cmath> /* for log2() */

#include "tree_distances.h"

using namespace Rcpp;

constexpr int32 LG2_SIZE = (int32(SL_MAX_TIPS) - 1) * (SL_MAX_TIPS - 1) + 1;

double lg2[LG2_SIZE];
double lg2_double_factorial[SL_MAX_TIPS + SL_MAX_TIPS - 2];
double lg2_unrooted[SL_MAX_TIPS + 2];
double *lg2_rooted = &lg2_unrooted[0] + 1;
__attribute__((constructor))
  void initialize_ldf() {
    lg2[0] = 0;
    for (int32 i = 1; i != LG2_SIZE; ++i) {
      lg2[i] = log2(i);
    }
    for (int16 i = 0; i != 3; ++i) {
      lg2_double_factorial[i] = 0;
      lg2_unrooted[i] = 0;
    }
    assert(lg2_rooted[0] == 0);
    assert(lg2_rooted[1] == 0);
    for (int32 i = 2; i != SL_MAX_TIPS + SL_MAX_TIPS - 2; ++i) {
      lg2_double_factorial[i] = lg2_double_factorial[i - 2] + lg2[i];
    }
    for (int32 i = 3; i != SL_MAX_TIPS + 2; ++i) {
      lg2_unrooted[i] = lg2_double_factorial[i + i - 5];
      assert(lg2_rooted[i - 1] == lg2_double_factorial[i + i - 5]);
    }
  }
