#ifndef _TREEDIST_TREE_DISTANCES_H
#define _TREEDIST_TREE_DISTANCES_H

#include "lap.h"
#include <TreeDist/mutual_clustering.h>
#include <TreeTools/SplitList.h>

/***** Re-export shared lookup tables to global scope *****/

using TreeDist::lg2;
using TreeDist::lg2_double_factorial;
using TreeDist::lg2_unrooted;
using TreeDist::lg2_rooted;

/***** Constants *****/

constexpr splitbit ALL_ONES = (std::numeric_limits<splitbit>::max)();

namespace TreeDist {

  // Re-exported from mutual_clustering.h:
  //   ic_matching(split_int a, split_int b, split_int n)

  void check_ntip(const int32 n);

  // See equation 16 in Meila 2007 (k' denoted K).
  // nkK is converted to pkK in the calling function when divided by n.
  inline void add_ic_element(double& ic_sum, const split_int nkK,
                             const split_int nk, const split_int nK,
                             const split_int n_tips,
                             const double lg2_n) noexcept {
    if (nkK && nk && nK) {
      ASSERT(n_tips > 0 && "n_tips must be positive");
      ASSERT(n_tips <= (std::numeric_limits<int32>::max)() / n_tips &&
        "nTip too large for int32 products in add_ic_element");
      assert(!(nkK == nk && nkK == nK && nkK << 1 == n_tips));
      const int32 numerator = nkK * n_tips;
      const int32 denominator = nk * nK;
      if (numerator != denominator) {
        ic_sum += nkK * (lg2_lookup(nkK) + lg2_n - lg2_lookup(nk) -
          lg2_lookup(nK));
      }
    }
  }


  [[nodiscard]] inline bool can_use_lookup_table(const split_int n_tips) noexcept {
    return n_tips <= static_cast<split_int>(SL_MAX_TIPS);
  }

  // Returns lg2_unrooted[x] - lg2_trees_matching_split(y, x - y)
  [[nodiscard]] inline double mmsi_pair_score_table(const split_int x,
                                                    const split_int y) noexcept {
    ASSERT(can_use_lookup_table(x));
    return lg2_unrooted[x] - (lg2_rooted[y] + lg2_rooted[x - y]);
  }

  [[nodiscard]] inline double mmsi_pair_score(const split_int x,
                                              const split_int y) noexcept {
    return lg2_unrooted_lookup(x) - (lg2_rooted_lookup(y) +
      lg2_rooted_lookup(x - y));
  }

  [[nodiscard]] inline double mmsi_score_table(const split_int n_same,
                                               const split_int n_a_and_b,
                                               const split_int n_different,
                                               const split_int n_a_only) noexcept {
    if (n_same == 0 || n_same == n_a_and_b)
      return mmsi_pair_score_table(n_different, n_a_only);
    if (n_different == 0 || n_different == n_a_only)
      return mmsi_pair_score_table(n_same, n_a_and_b);

    const double score1 = mmsi_pair_score_table(n_same, n_a_and_b),
      score2 = mmsi_pair_score_table(n_different, n_a_only);

    return (score1 > score2) ? score1 : score2;
  }

  [[nodiscard]] inline double mmsi_score(const split_int n_same,
                                         const split_int n_a_and_b,
                                         const split_int n_different,
                                         const split_int n_a_only) noexcept {
    if (n_same == 0 || n_same == n_a_and_b)
      return mmsi_pair_score(n_different, n_a_only);
    if (n_different == 0 || n_different == n_a_only)
      return mmsi_pair_score(n_same, n_a_and_b);

    const double score1 = mmsi_pair_score(n_same, n_a_and_b),
      score2 = mmsi_pair_score(n_different, n_a_only);

    return (score1 > score2) ? score1 : score2;
  }

  [[nodiscard]] inline double one_overlap_table(const split_int a, const split_int b,
                                                const split_int n) noexcept {
    ASSERT(can_use_lookup_table(n));
    if (a == b) {
      return lg2_rooted[a] + lg2_rooted[n - a];
    }
    const split_int lo = (a < b) ? a : b;
    const split_int hi = (a < b) ? b : a;
    return lg2_rooted[hi] + lg2_rooted[n - lo] - lg2_rooted[hi - lo + 1];
  }

  [[nodiscard]] inline double one_overlap(const split_int a, const split_int b,
                                          const split_int n) noexcept {
    if (a == b) {
      return lg2_rooted_lookup(a) + lg2_rooted_lookup(n - a);
    }
    const split_int lo = (a < b) ? a : b;
    const split_int hi = (a < b) ? b : a;
    return lg2_rooted_lookup(hi) + lg2_rooted_lookup(n - lo) -
      lg2_rooted_lookup(hi - lo + 1);
  }

  [[nodiscard]] inline double one_overlap_notb_table(const split_int a,
                                                      const split_int n_minus_b,
                                                      const split_int n) noexcept {
    ASSERT(can_use_lookup_table(n));
    const split_int b = n - n_minus_b;
    if (a == b) {
      return lg2_rooted[b] + lg2_rooted[n_minus_b];
    } else if (a < b) {
      return lg2_rooted[b] + lg2_rooted[n - a] - lg2_rooted[b - a + 1];
    } else {
      return lg2_rooted[a] + lg2_rooted[n_minus_b] - lg2_rooted[a - b + 1];
    }
  }

  [[nodiscard]] inline double one_overlap_notb(const split_int a,
                                               const split_int n_minus_b,
                                               const split_int n) noexcept {
    const split_int b = n - n_minus_b;
    if (a == b) {
      return lg2_rooted_lookup(b) + lg2_rooted_lookup(n_minus_b);
    } else if (a < b) {
      return lg2_rooted_lookup(b) + lg2_rooted_lookup(n - a) -
        lg2_rooted_lookup(b - a + 1);
    } else {
      return lg2_rooted_lookup(a) + lg2_rooted_lookup(n_minus_b) -
        lg2_rooted_lookup(a - b + 1);
    }
  }

  [[nodiscard]] inline double spi_overlap_table(const splitbit* a_state,
                                                const splitbit* b_state,
                                                const split_int n_tips,
                                                const split_int in_a,
                                                const split_int in_b,
                                                const split_int n_bins) noexcept {
    ASSERT(can_use_lookup_table(n_tips));
    split_int n_ab = 0;
    for (split_int bin = 0; bin < n_bins; ++bin) {
      n_ab += TreeTools::count_bits(a_state[bin] & b_state[bin]);
    }

    if (n_ab == 0) {
      return one_overlap_notb_table(in_a, in_b, n_tips);
    }
    if (n_ab == in_b || n_ab == in_a) {
      return one_overlap_table(in_a, in_b, n_tips);
    }
    if (in_a + in_b - n_ab == n_tips) {
      return one_overlap_notb_table(in_a, in_b, n_tips);
    }

    return 0.0;
  }

  // Popcount-based: single pass over bins replaces 4 sequential boolean scans.
  // Counts n_ab = |A ∩ B| via hardware POPCNT, then derives all 4 Venn-diagram
  // region populations from arithmetic on n_ab, in_a, in_b, n_tips.
  [[nodiscard]] inline double spi_overlap(const splitbit* a_state, const splitbit* b_state,
                                          const split_int n_tips, const split_int in_a,
                                          const split_int in_b, const split_int n_bins) noexcept {
    split_int n_ab = 0;
    for (split_int bin = 0; bin < n_bins; ++bin) {
      n_ab += TreeTools::count_bits(a_state[bin] & b_state[bin]);
    }

    if (n_ab == 0) {
      return one_overlap_notb(in_a, in_b, n_tips);
    }
    if (n_ab == in_b || n_ab == in_a) {
      return one_overlap(in_a, in_b, n_tips);
    }
    if (in_a + in_b - n_ab == n_tips) {
      return one_overlap_notb(in_a, in_b, n_tips);
    }

    return 0.0;
  }
}

extern Rcpp::List cpp_robinson_foulds_distance(Rcpp::RawMatrix x,
                                               Rcpp::RawMatrix y,
                                               Rcpp::IntegerVector nTip);

#endif
