#include <Rcpp/Lightest>
#include <BigTreeTools/SplitList.h> /* for SL_MAX_TIPS */

#include <cmath> /* for log2() */

#include "tree_distances.h"

using namespace Rcpp;

#define LG2_SIZE ((int32(SL_MAX_TIPS) - 1) * (SL_MAX_TIPS - 1) + 1)

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


// Returns lg2_unrooted[x] - lg2_trees_matching_split(y, x - y)
double mmsi_pair_score (const int16 x, const int16 y) {
  assert(SL_MAX_TIPS + 2 <= INT_16_MAX); // verify int16 ok
  
  return lg2_unrooted[x] - (lg2_rooted[y] + lg2_rooted[x - y]);
}

double mmsi_score(const int16 n_same, const int16 n_a_and_b,
                  const int16 n_different, const int16 n_a_only) {
  if (n_same == 0 || n_same == n_a_and_b)
    return mmsi_pair_score(n_different, n_a_only);
  if (n_different == 0 || n_different == n_a_only)
    return mmsi_pair_score(n_same, n_a_and_b);
  
  const double
    score1 = mmsi_pair_score(n_same, n_a_and_b),
      score2 = mmsi_pair_score(n_different, n_a_only);
  
  return (score1 > score2) ? score1 : score2;
}


/* 
 * See equation 16 in Meila 2007. I denote k' as K.
 * nkK is converted to pkK in the calling function, when the sum of all
 * elements is divided by n.
*/
double ic_element (const int16 nkK, const int16 nk,
                   const int16 nK, const int16 n) {
  if (nkK && nk && nK) {
    // Avoid possible rounding errors
    
    if (nkK == nk && nkK == nK && nkK + nkK == n) return nkK;
    const int32
      numerator = nkK * n,
      denominator = nk * nK
    ;
    if (numerator == denominator) return 0; 
    
    // Multiply-then-log is twice as fast log-then-add
    return nkK * (lg2[numerator] - lg2[denominator]);
  } else return 0;
}

double ic_matching (const int16 a, const int16 b, const int16 n) {
  return (a * (lg2[n] - lg2[a])) + 
         (b * (lg2[n] - lg2[b]));
}

double one_overlap (const int16 a, const int16 b, const int16 n) {
  assert(SL_MAX_TIPS + 2 <= INT_16_MAX); // verify int16 ok
  if (a == b) {
    return lg2_rooted[a] + lg2_rooted[n - a];
  } else if (a < b) {
    return lg2_rooted[b] + lg2_rooted[n - a] - lg2_rooted[b - a + 1];
  } else {
    return lg2_rooted[a] + lg2_rooted[n - b] - lg2_rooted[a - b + 1];
  }
}

double one_overlap_notb (const int16 a, const int16 n_minus_b, const int16 n) {
  assert(SL_MAX_TIPS + 2 <= INT_16_MAX); // verify int16 ok
  const int16 b = n - n_minus_b;
  if (a == b) {
    return lg2_rooted[b] + lg2_rooted[n_minus_b];
  } else if (a < b) {
    return lg2_rooted[b] + lg2_rooted[n - a] - lg2_rooted[b - a + 1];
  } else {
    return lg2_rooted[a] + lg2_rooted[n_minus_b] - lg2_rooted[a - b + 1];
  }
}

double spi_overlap (const splitbit* a_state, const splitbit* b_state,
                    const int16 n_tips, const int16 in_a, const int16 in_b,
                    const int16 n_bins) {
  bool flag = true;
  
  assert(SL_MAX_BINS <= INT16_MAX);
  for (int16 bin = 0; bin != n_bins; bin++) {
    if (a_state[bin] & b_state[bin]) {
      flag = false;
      break;
    }
  }
  if (flag) return one_overlap_notb(in_a, in_b, n_tips);
  
  for (int16 bin = 0; bin != n_bins; bin++) {
    if ((~a_state[bin] & b_state[bin])) {
      flag = true;
      break;
    }
  }
  if (!flag) return one_overlap(in_a, in_b, n_tips);
  
  for (int16 bin = 0; bin != n_bins; bin++) {
    if ((a_state[bin] & ~b_state[bin])) {
      flag = false;
      break;
    }
  }
  if (flag) return one_overlap(in_a, in_b, n_tips);
  
  const int16 loose_end_tips = n_tips % SL_BIN_SIZE;
  for (int16 bin = 0; bin != n_bins; bin++) {
    splitbit test = ~(a_state[bin] | b_state[bin]);
    if (bin == n_bins - 1 && loose_end_tips) {
      test &= ~(ALL_ONES << loose_end_tips);
    }
    if (test) {
      flag = true;
      break;
    }
  }
  if (!flag) return one_overlap_notb(in_a, in_b, n_tips);
  
  return 0;
}
