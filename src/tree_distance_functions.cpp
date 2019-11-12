#include <Rcpp.h>
using namespace Rcpp;
#include <stdint.h>
#include <math.h> /* for pow(), log() */
#include "tree_distances.h"
#include "SplitList.h"


uint32_t bitcounts[65536]; // the bytes representing bit count of each number 0-65535
__attribute__((constructor))
  void initialize_bitcounts() {
    for (int32_t i = 0; i < 65536; i++) {
      int32_t n_bits = 0;
      for (int j = 0; j < 16; j++) {
        if ((i & powers_of_two[j])) ++n_bits;
      }
      bitcounts[i] = n_bits;
    }
  }

int count_bits (splitbit x) {
  /* For 32-bit splitbits: */
  /* return bitcounts[x & right16bits] + bitcounts[x >> 16]; */
  
  /* For 64-bit splitbits: */
  return bitcounts[x & right16bits] + bitcounts[(x >> 16) & right16bits]
  + bitcounts[(x >> 32) & right16bits]
  + bitcounts[(x >> 48)];
}

double lg2_double_factorial[MAX_TIPS + MAX_TIPS - 2];
double lg2_rooted[MAX_TIPS + 1];
double lg2_unrooted[MAX_TIPS + 1];
__attribute__((constructor))
  void initialize_ldf() {
    for (int i = 0; i != 3; i++) {
      lg2_double_factorial[i] = 0;
      lg2_rooted[i] = 0;
      lg2_unrooted[i] = 0;
    }
    for (int i = 2; i != MAX_TIPS + MAX_TIPS - 2; i++) {
      lg2_double_factorial[i] = lg2_double_factorial[i - 2] + log2(i);
    }
    for (int i = 3; i != MAX_TIPS + 1; i++) {
      lg2_unrooted[i] = lg2_double_factorial[i + i - 5];
      lg2_rooted[i] = lg2_double_factorial[i + i - 3];
    }
  }

double lg2_trees_matching_split (int a, int b) {
  if (a == 0) return lg2_unrooted[b];
  if (b == 0) return lg2_unrooted[a];
  return lg2_rooted[a] + lg2_rooted[b];
}


double p_lg2_p_frac (double p) {
  return -p * log2(p);
}

double p_lg2_p (double p) {
  if (p == 0) return 0;
  if (p == 1) return 0;
  return p_lg2_p_frac(p);
}

double entropy2 (double p) {
  if (p == 0) return 0;
  if (p == 1) return 0;
  return p_lg2_p_frac(p) + p_lg2_p_frac(1 - p);
}

double entropy4 (double p1, double p2, double p3, double p4) {
  return p_lg2_p(p1) +  p_lg2_p(p2) +  p_lg2_p(p3) +  p_lg2_p(p4);
}


double one_overlap (const int *a, const int *b, const int *n) {
  if (*a == *b) return lg2_rooted[*a] + lg2_rooted[*n - *a];
  if (*a < *b) return lg2_rooted[*b] + lg2_rooted[*n -*a] - lg2_rooted[*b - *a + 1];
  return lg2_rooted[*a] + lg2_rooted[*n - *b] - lg2_rooted[*a - *b + 1];
}

double one_overlap_notb (const int *a, const int *n_minus_b, const int *n) {
  const int b = *n - *n_minus_b;
  if (*a == b) return lg2_rooted[b] + lg2_rooted[*n_minus_b];
  if (*a < b) return lg2_rooted[b] + lg2_rooted[*n - *a] - lg2_rooted[b -* a + 1];
  return lg2_rooted[*a] + lg2_rooted[*n_minus_b] - lg2_rooted[*a - b + 1];
}

double mpi (const splitbit* a_state, const splitbit* b_state,
            const int *n_tips, const int *in_a, const int *in_b,
            const double *lg2_unrooted_n, const int *n_bins) {
  bool flag = true;
  
    for (int bin = 0; bin < *n_bins; bin++) {
    if (a_state[bin] & b_state[bin]) {
      flag = false;
      break;
    }
  }
  if (flag) return *lg2_unrooted_n - one_overlap_notb(in_a, in_b, n_tips);
  
  for (int bin = 0; bin < *n_bins; bin++) {
    if ((~a_state[bin] & b_state[bin])) {
      flag = true;
      break;
    }
  }
  if (!flag) return *lg2_unrooted_n - one_overlap(in_a, in_b, n_tips);
  
  for (int bin = 0; bin < *n_bins; bin++) {
    if ((a_state[bin] & ~b_state[bin])) {
      flag = false;
      break;
    }
  }
  if (flag) return *lg2_unrooted_n - one_overlap(in_a, in_b, n_tips);
  
  for (int bin = 0; bin < *n_bins; bin++) {
    splitbit test = ~(a_state[bin] | b_state[bin]);
    if (bin == *n_bins - 1) test &= ~(ALL_ONES << (*n_tips % BIN_SIZE));
    if (test) {
      flag = true;
      break;
    }
  }
  if (!flag) return *lg2_unrooted_n - one_overlap_notb(in_a, in_b, n_tips);
  
  return 0;
}
