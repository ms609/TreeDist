#ifndef TREEDIST_MUTUAL_CLUSTERING_H_
#define TREEDIST_MUTUAL_CLUSTERING_H_

// Mutual Clustering Information (MCI) — public API.
//
// Provides the core MCI score computation between two split sets,
// clustering entropy, and lookup-table initialization.
//
// Implementation lives in mutual_clustering_impl.h, guarded by
// TREEDIST_MCI_IMPLEMENTATION.

#include "lap.h"
#include <TreeTools/SplitList.h>   // splitbit, SL_MAX_TIPS, count_bits
#include <cmath>

namespace TreeDist {

  // ---- Lookup tables (populated by init_lg2_tables) ----
  //
  // lg2[i]                = log2(i)  for  0 <= i <= SL_MAX_TIPS
  // lg2_double_factorial[i] = log2(i!!) for 0 <= i < 2*SL_MAX_TIPS - 2
  // lg2_unrooted[i]       = log2((2i-5)!!) for i >= 3
  // lg2_rooted             = &lg2_unrooted[0] + 1  (so lg2_rooted[i] = lg2_unrooted[i+1])
  //
  // These are fast-path caches sized to TreeTools' stack threshold.
  // For larger trees, fall back to on-the-fly computation via the _lookup helpers.
  // Definitions are in mutual_clustering_impl.h.

  extern double lg2[SL_MAX_TIPS + 1];
  extern double lg2_double_factorial[SL_MAX_TIPS + SL_MAX_TIPS - 2];
  extern double lg2_unrooted[SL_MAX_TIPS + 2];
  extern double* lg2_rooted;

  // Populate lookup tables.  Must be called once before any MCI
  // computation.  max_tips should be >= the largest tree size used.
  void init_lg2_tables(int max_tips);

  // Out-of-line slow paths for inputs that exceed the precomputed tables
  // (i.e. n_tips > SL_MAX_TIPS).  Defined in mutual_clustering_impl.h.
  // Kept out-of-line so the *_lookup() inline thunks below stay tiny —
  // table load + predicted branch, no extern function calls — which some
  // compilers refuse to inline through, suppressing optimization in the
  // per-cell hot loops of PID/MCI/MSI/RF info.
  double lg2_slow(split_int x);
  double lg2_unrooted_slow(split_int n_tips);
  double lg2_rooted_slow(split_int n_tips);

  // log2(x) — table-fast for x <= SL_MAX_TIPS, runtime otherwise.
  [[nodiscard]] inline double lg2_lookup(split_int x) noexcept {
    return (x <= static_cast<split_int>(SL_MAX_TIPS))
      ? lg2[x]
      : lg2_slow(x); // LCOV_EXCL_LINE
  }

  // log2((2n-5)!!) — table-fast for n <= SL_MAX_TIPS+1, lgamma otherwise.
  [[nodiscard]] inline double lg2_unrooted_lookup(split_int n_tips) noexcept {
    return (n_tips <= static_cast<split_int>(SL_MAX_TIPS + 1))
      ? lg2_unrooted[n_tips]
      : lg2_unrooted_slow(n_tips); // LCOV_EXCL_LINE
  }

  // log2((2n-3)!!) — table-fast for n <= SL_MAX_TIPS+1, lgamma otherwise.
  [[nodiscard]] inline double lg2_rooted_lookup(split_int n_tips) noexcept {
    return (n_tips <= static_cast<split_int>(SL_MAX_TIPS + 1))
      ? lg2_rooted[n_tips]
      : lg2_rooted_slow(n_tips); // LCOV_EXCL_LINE
  }

  // ---- Inline helpers ----

  // Information content of a perfectly-matching split pair.
  // ic_matching(a, b, n) = (a + b) * lg2[n] - a * lg2[a] - b * lg2[b]
  [[nodiscard]] inline double ic_matching(split_int a, split_int b,
                                          split_int n) noexcept {
    const double lg2a = lg2_lookup(a);
    const double lg2b = lg2_lookup(b);
    const double lg2n = lg2_lookup(n);
    return (a + b) * lg2n - a * lg2a - b * lg2b;
  }

  // Clustering entropy of a split set.
  // Computes H_clust = -sum_i [ p_i * log2(p_i) + q_i * log2(q_i) ]
  // where p_i = in_split[i] / n_tips and q_i = 1 - p_i.
  [[nodiscard]] inline double clustering_entropy(
      const int* in_split, int n_splits, int n_tips) {
    if (n_tips <= 1 || n_splits == 0) return 0.0;
    const double lg2n = std::log2(static_cast<double>(n_tips));
    double ce = 0.0;
    for (int i = 0; i < n_splits; ++i) {
      const int a = in_split[i];
      const int b = n_tips - a;
      if (a > 1 && b > 1) {
        ce += (a * std::log2(static_cast<double>(a)) +
               b * std::log2(static_cast<double>(b)))
              / n_tips - lg2n;
      }
    }
    return -ce;
  }

  // ---- MCI score (declaration only) ----
  //
  // Computes the Mutual Clustering Information between two split sets.
  // Uses sort+merge exact-match detection + LAP on the reduced matrix.
  //
  // a_state[i] points to n_bins splitbit words for split i of tree A.
  // a_in[i]    = popcount of split i (tips in the "1" partition).
  // Returns the MCI score (higher = more similar).
  //
  // Implementation in mutual_clustering_impl.h.

  double mutual_clustering_score(
      const splitbit* const* a_state, const split_int* a_in, split_int a_n_splits,
      const splitbit* const* b_state, const split_int* b_in, split_int b_n_splits,
      split_int n_bins, int32 n_tips,
      LapScratch& scratch);

} // namespace TreeDist

#endif // TREEDIST_MUTUAL_CLUSTERING_H_
