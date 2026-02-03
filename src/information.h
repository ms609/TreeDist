#ifndef _TREEDIST_INFO_H
#define _TREEDIST_INFO_H

#include <cmath> /* for log2() */
#include <TreeTools/ClusterTable.h> /* for CT_MAX_LEAVES */

#include "ints.h" /* for int16 */

constexpr int_fast32_t LOG_MAX = 2048;
static double log2_table[LOG_MAX];
__attribute__((constructor))
void compute_log2_table() {
  for (int i = 1; i < LOG_MAX; ++i) {
    log2_table[i] = std::log2(i);
  }
}

constexpr int_fast32_t FACT_MAX = CT_MAX_LEAVES + CT_MAX_LEAVES + 5 + 1;
constexpr double log_2 = 0.6931471805599452862268;

double ldfact[FACT_MAX];
double l2rooted[CT_MAX_LEAVES];
double *l2unrooted = &l2rooted[0] - 1;
double l2[CT_MAX_LEAVES];
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
    
    for (int_fast32_t i = 3; i != CT_MAX_LEAVES; ++i) {
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


//' Calculate entropy of integer vector of counts
//' 
//' Wrapper for C++ function; no input checking is performed.
//' [`Ntropy()`] is better suited for use where performance is not critical.
//' @param n a vector of integer counts
//' @return `entropy_int()` returns a numeric corresponding to the entropy of
//' each observation, in bits.
//' @export
//' @keywords internal
// [[Rcpp::export]]
double entropy_int(const Rcpp::IntegerVector &n) {
  int N = 0;
  double sum = 0;
  
  for (int x : n) {
    if (x > 0) {
      if (x < LOG_MAX) {
        sum += x * log2_table[x];
      } else {
        sum += x * std::log2(x);
      }
      N += x;
    }
  }
  
  return (N == 0) ? 0.0 : (N < LOG_MAX ? log2_table[N] : std::log2(N)) - sum / N;
}


// TODO This is copied from TreeSearch/src/expected_mi.cpp - double definition
// is bad practice; define once (here? TreeTools?) and cross-reference.
// TODO this also somewhat duplicates the above...
#define MAX_FACTORIAL_LOOKUP 8192
static double log2_factorial_table[MAX_FACTORIAL_LOOKUP + 1];
static const double LOG2_E = 1.4426950408889634;

__attribute__((constructor))
void initialize_factorial_cache() {
  log2_factorial_table[0] = 0.0;
  for (int i = 1; i <= MAX_FACTORIAL_LOOKUP; i++) {
    log2_factorial_table[i] = log2_factorial_table[i - 1] + std::log2(i);
  }
}
// Fast lookup with bounds checking
inline double l2factorial(int n) {
  if (n <= MAX_FACTORIAL_LOOKUP) {
    return log2_factorial_table[n];
  } else {
    return lgamma(n + 1) * LOG2_E;
  }
}

// ni and nj are vectors listing the number of entitites in each cluster
// [[Rcpp::export]]
double expected_mi(const IntegerVector &ni, const IntegerVector &nj) {
  // ni = {a, N-a}; nj = counts of character states
  const int a = ni[0];
  const int N = ni[0] + ni[1];
  if (a <= 0 || a >= N) return 0.0; // trivial split
  
  const double invN = 1.0 / static_cast<double>(N);
  const double log2N = std::log2(static_cast<double>(N));
  const double log2a  = std::log2(static_cast<double>(a));
  const double log2Na = std::log2(static_cast<double>(N - a));
  const double log2_denom = l2factorial(N) - l2factorial(a) - l2factorial(N - a);
  
  double emi = 0.0;
  
  for (int j = 0; j < nj.size(); ++j) {
    int mj = nj[j];
    if (mj <= 0) continue;
    
    int kmin = std::max(0, a + mj - N);
    int kmax = std::min(a, mj);
    if (kmin > kmax) continue;
    
    const double log2mj = std::log2(static_cast<double>(mj));
    
    // compute P(K=kmin)
    double log2P = (l2factorial(mj) - l2factorial(kmin) - l2factorial(mj - kmin))
      + (l2factorial(N - mj) - l2factorial(a - kmin) - l2factorial(N - mj - (a - kmin)))
      - log2_denom;
      double Pk = std::pow(2.0, log2P);
      
      for (int k = kmin; k <= kmax; ++k) {
        if (Pk > 0.0) {
          // contribution from inside the split
          if (k > 0) {
            double mi_in = std::log2(static_cast<double>(k)) + log2N - (log2a + log2mj);
            emi += (static_cast<double>(k) * invN) * mi_in * Pk;
          }
          // contribution from outside the split
          int kout = mj - k;
          if (kout > 0) {
            double mi_out = std::log2(static_cast<double>(kout)) + log2N - (log2Na + log2mj);
            emi += (static_cast<double>(kout) * invN) * mi_out * Pk;
          }
        }
        // Update P(k) → P(k+1)
        if (k < kmax) {
          double numer = static_cast<double>((mj - k) * (a - k));
          double denom = static_cast<double>((k + 1) * (N - mj - a + k + 1));
          Pk *= numer / denom;
        }
      }
  }
  
  return emi;
}


#endif

