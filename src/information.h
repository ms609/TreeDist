#ifndef _TREEDIST_INFO_H
#define _TREEDIST_INFO_H

#include <cmath> /* for log2() */
#include "ints.h" /* for int16 */

const int_fast32_t
  DAY_MAX_LEAVES = 16383,
  FACT_MAX = DAY_MAX_LEAVES + DAY_MAX_LEAVES + 5 + 1
;

const double log_2 = log(2.0);
double ldfact[FACT_MAX];
double l2rooted[DAY_MAX_LEAVES];
double *l2unrooted = &l2rooted[0] - 1;
double l2[DAY_MAX_LEAVES];
__attribute__((constructor))
  void compute_double_factorials() {
    ldfact[0] = 0;
    ldfact[1] = 0;
    l2rooted[0] = 0;
    l2rooted[1] = 0;
    l2rooted[2] = 0;
    l2[1] = 0;
    l2[2] = 1;
    
    for (int_fast32_t i = 2; i != FACT_MAX; ++i) {
      ldfact[i] = ldfact[i - 2] + log2(double(i));
    }
    
    for (int_fast32_t i = 3; i != DAY_MAX_LEAVES; ++i) {
      l2[i] = log2(double(i));
      l2rooted[i] = ldfact[(i << 1) - 3];
      assert(l2unrooted[i] == ldfact[(i << 1) - 5]);
    }
  }

inline double split_phylo_info (const int16 n_in, const int16 *n_tip,
                                const double p) {
  const int16 n_out = *n_tip - n_in;
  assert(p > 0);
  assert(p <= 1);
  assert(n_in > 1);
  assert(n_out > 1);
  if (p == 1) {
    return (l2unrooted[*n_tip] - l2rooted[n_in] - l2rooted[n_out]);
  } else {
    const double 
      q = 1 - p,
      l2n = l2unrooted[*n_tip],
      l2n_consistent = l2rooted[n_in] + l2rooted[n_out],
      l2p_consistent = l2n_consistent - l2n,
      l2p_inconsistent = log2(-expm1(l2p_consistent * log_2)),
      l2n_inconsistent = l2p_inconsistent + l2n
    ;
    
    return(l2n +
           p * (log2(p) - l2n_consistent) +
           q * (log2(q) - l2n_inconsistent));
  }
}

inline double split_clust_info (const int16 n_in, const int16 *n_tip,
                                const double p) {
  const int16 n_out = *n_tip - n_in;
  assert(p > 0);
  assert(p <= 1);
  assert(n_in > 1);
  assert(n_out > 1);
  return -p * (
      (((n_in * l2[n_in]) + (n_out * l2[n_out])) / *n_tip) - l2[*n_tip]
  );
}

#endif