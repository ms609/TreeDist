#ifndef _TREEDIST_TREE_DISTANCES_H
#define _TREEDIST_TREE_DISTANCES_H

#include "lap.h"
#include <TreeDist/mutual_clustering.h>
#include <TreeTools/SplitList.h>

/***** Re-export shared lookup tables and fallback functions to global scope *****/

using TreeDist::lg2;
using TreeDist::lg2_double_factorial;
using TreeDist::lg2_unrooted;
using TreeDist::lg2_rooted;
using TreeDist::lg2_lookup;
using TreeDist::lg2_unrooted_lookup;
using TreeDist::lg2_rooted_lookup;

/***** Constants *****/

constexpr splitbit ALL_ONES = (std::numeric_limits<splitbit>::max)();

namespace TreeDist {

  // Validate that n_tips does not exceed the compiled SL_MAX_TIPS limit.
  // Defined in tree_distances.cpp; calls Rcpp::stop() on failure.
  void check_ntip(int32 n);

  // Re-exported from mutual_clustering.h:
  //   ic_matching(split_int a, split_int b, split_int n)

  // See equation 16 in Meila 2007 (k' denoted K).
  // nkK is converted to pkK in the calling function when divided by n.
  // Parameters are split_int (int32) to handle n_tips up to 32767 safely;
  // products nkK*n_tips and nk*nK fit int32 since both factors <= 32767.
  inline void add_ic_element(double& ic_sum, const split_int nkK, const split_int nk,
                             const split_int nK, const split_int n_tips,
                             const double lg2_n) noexcept {
    if (nkK && nk && nK) {
      assert(!(nkK == nk && nkK == nK && nkK << 1 == n_tips));
      const int32 numerator = nkK * n_tips;
      const int32 denominator = nk * nK;
      if (numerator != denominator) {
        ic_sum += nkK * (lg2_lookup(nkK) + lg2_n - lg2_lookup(nk) - lg2_lookup(nK));
      }
    }
  }


  // Returns lg2_unrooted[x] - lg2_trees_matching_split(y, x - y)
  [[nodiscard]] inline double mmsi_pair_score(const split_int x, const split_int y) noexcept {
    static_assert(SL_MAX_TIPS + 2 <= std::numeric_limits<int16>::max(),
                  "int16 too narrow for SL_MAX_TIPS");

    return lg2_unrooted_lookup(x) - (lg2_rooted_lookup(y) + lg2_rooted_lookup(x - y));
  }

  [[nodiscard]] inline double mmsi_score(const split_int n_same, const split_int n_a_and_b,
                    const split_int n_different, const split_int n_a_only)  noexcept {
    if (n_same == 0 || n_same == n_a_and_b)
      return mmsi_pair_score(n_different, n_a_only);
    if (n_different == 0 || n_different == n_a_only)
      return mmsi_pair_score(n_same, n_a_and_b);

    const double
      score1 = mmsi_pair_score(n_same, n_a_and_b),
        score2 = mmsi_pair_score(n_different, n_a_only);

    return (score1 > score2) ? score1 : score2;
  }


  // Compile-time-dispatched lg2 helpers: when Fast=true, callers guarantee
  // x <= SL_MAX_TIPS, so the bounds-check + branch in lg2_rooted_lookup is
  // skipped entirely and the inner kernel collapses to a bare table load.
  // Used by spi_overlap / one_overlap below; large-tree callers pass
  // Fast=false (default) and pay the predictable branch.
  template<bool Fast>
  [[nodiscard]] inline double lg2_rooted_fast(split_int x) noexcept {
    if constexpr (Fast) return lg2_rooted[x];
    else return lg2_rooted_lookup(x);
  }

  template<bool Fast>
  [[nodiscard]] inline double lg2_unrooted_fast(split_int x) noexcept {
    if constexpr (Fast) return lg2_unrooted[x];
    else return lg2_unrooted_lookup(x);
  }

  template<bool Fast = false>
  [[nodiscard]] inline double one_overlap(const split_int a, const split_int b, const split_int n) noexcept {
    static_assert(SL_MAX_TIPS + 2 <= std::numeric_limits<int16>::max(),
                  "int16 too narrow for SL_MAX_TIPS");
    if (a == b) {
      return lg2_rooted_fast<Fast>(a) + lg2_rooted_fast<Fast>(n - a);
    }
    // Unify a<b and a>b via lo/hi: removes an unpredictable branch.
    const split_int lo = (a < b) ? a : b;
    const split_int hi = (a < b) ? b : a;
    return lg2_rooted_fast<Fast>(hi) + lg2_rooted_fast<Fast>(n - lo)
         - lg2_rooted_fast<Fast>(hi - lo + 1);
  }

  template<bool Fast = false>
  [[nodiscard]] inline double one_overlap_notb(const split_int a, const split_int n_minus_b, const split_int n) noexcept {
    static_assert(SL_MAX_TIPS + 2 <= std::numeric_limits<int16>::max(),
                  "int16 too narrow for SL_MAX_TIPS");
    const split_int b = n - n_minus_b;
    if (a == b) {
      return lg2_rooted_fast<Fast>(b) + lg2_rooted_fast<Fast>(n_minus_b);
    } else if (a < b) {
      return lg2_rooted_fast<Fast>(b) + lg2_rooted_fast<Fast>(n - a)
           - lg2_rooted_fast<Fast>(b - a + 1);
    } else {
      return lg2_rooted_fast<Fast>(a) + lg2_rooted_fast<Fast>(n_minus_b)
           - lg2_rooted_fast<Fast>(a - b + 1);
    }
  }


  // Popcount-based: single pass over bins replaces 4 sequential boolean scans.
  // Counts n_ab = |A ∩ B| via hardware POPCNT, then derives all 4 Venn-diagram
  // region populations from arithmetic on n_ab, in_a, in_b, n_tips.
  // split_int used throughout so in_a + in_b does not overflow (max ~65534 for 32767-tip trees).
  //
  // Fast=true: caller guarantees n_tips <= SL_MAX_TIPS, so inner one_overlap*
  // calls use direct lg2_rooted[] access.  Caller (shared_phylo_score) dispatches
  // once on n_tips; per-cell branches are eliminated.
  template<bool Fast = false>
  [[nodiscard]] inline double spi_overlap(const splitbit* a_state, const splitbit* b_state,
                       const split_int n_tips, const split_int in_a,
                       const split_int in_b, const split_int n_bins) noexcept {

    static_assert(SL_MAX_BINS <= INT16_MAX,
                  "int16 too narrow for SL_MAX_BINS");

    split_int n_ab = 0;
    for (split_int bin = 0; bin < n_bins; ++bin) {
      n_ab += TreeTools::count_bits(a_state[bin] & b_state[bin]);
    }

    // n_a_only  = in_a - n_ab   (tips in A but not B)
    // n_b_only  = in_b - n_ab   (tips in B but not A)
    // n_neither = n_tips - in_a - in_b + n_ab  (tips in neither)
    //
    // Return 0 when all 4 regions are populated (the common case for
    // unrelated splits).  Otherwise return the appropriate one_overlap score.

    if (n_ab == 0) {
      return one_overlap_notb<Fast>(in_a, in_b, n_tips);
    }
    if (n_ab == in_b || n_ab == in_a) {
      // B ⊆ A (n_b_only == 0) or A ⊆ B (n_a_only == 0)
      return one_overlap<Fast>(in_a, in_b, n_tips);
    }
    if (in_a + in_b - n_ab == n_tips) {
      // A ∪ B covers all tips (n_neither == 0)
      return one_overlap_notb<Fast>(in_a, in_b, n_tips);
    }

    return 0.0;
  }
}

extern Rcpp::List cpp_robinson_foulds_distance(Rcpp::RawMatrix x,
                                               Rcpp::RawMatrix y,
                                               Rcpp::IntegerVector nTip);

#endif
