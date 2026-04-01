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
  // These are defined in mutual_clustering_impl.h.

  extern double lg2[SL_MAX_TIPS + 1];
  extern double lg2_double_factorial[SL_MAX_TIPS + SL_MAX_TIPS - 2];
  extern double lg2_unrooted[SL_MAX_TIPS + 2];
  extern double* lg2_rooted;

  // Populate lookup tables.  Must be called once before any MCI
  // computation.  max_tips should be >= the largest tree size used.
  void init_lg2_tables(int max_tips);

  // ---- Inline helpers ----

  // Information content of a perfectly-matching split pair.
  // ic_matching(a, b, n) = (a + b) * lg2[n] - a * lg2[a] - b * lg2[b]
  [[nodiscard]] inline double ic_matching(int16 a, int16 b,
                                          int16 n) noexcept {
    const double lg2a = lg2[a];
    const double lg2b = lg2[b];
    const double lg2n = lg2[n];
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
      const splitbit* const* a_state, const int16* a_in, int16 a_n_splits,
      const splitbit* const* b_state, const int16* b_in, int16 b_n_splits,
      int16 n_bins, int32 n_tips,
      LapScratch& scratch);

} // namespace TreeDist

#endif // TREEDIST_MUTUAL_CLUSTERING_H_
