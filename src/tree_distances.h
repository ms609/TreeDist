#ifndef _TREEDIST_TREE_DISTANCES_H
#define _TREEDIST_TREE_DISTANCES_H

#include "lap.h"

/***** Constants requiring initialization *****/

constexpr splitbit ALL_ONES = (std::numeric_limits<splitbit>::max)();
extern double lg2[SL_MAX_TIPS + 1];
extern double lg2_double_factorial[SL_MAX_TIPS + SL_MAX_TIPS - 2];
extern double lg2_unrooted[SL_MAX_TIPS + 2];
extern double *lg2_rooted;

namespace TreeDist {

  // See equation 16 in Meila 2007 (k' denoted K).
  // nkK is converted to pkK in the calling function when divided by n.
  inline void add_ic_element(double& ic_sum, const int16 nkK, const int16 nk,
                             const int16 nK, const int16 n_tips,
                             const double lg2_n) noexcept {
    if (nkK && nk && nK) {
      assert(!(nkK == nk && nkK == nK && nkK << 1 == n_tips));
      const int32 numerator = nkK * n_tips;
      const int32 denominator = nk * nK;
      if (numerator != denominator) {
        ic_sum += nkK * (lg2[nkK] + lg2_n - lg2[nk] - lg2[nK]);
      }
    }
  }


  // Returns lg2_unrooted[x] - lg2_trees_matching_split(y, x - y)
  [[nodiscard]] inline double mmsi_pair_score(const int16 x, const int16 y) noexcept {
    assert(SL_MAX_TIPS + 2 <= std::numeric_limits<int16>::max()); // verify int16 ok
    
    return lg2_unrooted[x] - (lg2_rooted[y] + lg2_rooted[x - y]);
  }

  [[nodiscard]] inline double mmsi_score(const int16 n_same, const int16 n_a_and_b,
                    const int16 n_different, const int16 n_a_only)  noexcept {
    if (n_same == 0 || n_same == n_a_and_b)
      return mmsi_pair_score(n_different, n_a_only);
    if (n_different == 0 || n_different == n_a_only)
      return mmsi_pair_score(n_same, n_a_and_b);
    
    const double
      score1 = mmsi_pair_score(n_same, n_a_and_b),
        score2 = mmsi_pair_score(n_different, n_a_only);
    
    return (score1 > score2) ? score1 : score2;
  }


[[nodiscard]] inline double ic_matching(const int16 a, const int16 b, const int16 n) noexcept {
    const double lg2a = lg2[a];
    const double lg2b = lg2[b];
    const double lg2n = lg2[n];
    return (a + b) * lg2n - a * lg2a - b * lg2b;
    //  (a * (lg2n - lg2a)) + (b * (lg2n - lg2b)); is substantially slower
  }

[[nodiscard]] inline double one_overlap(const int16 a, const int16 b, const int16 n) noexcept {
    assert(SL_MAX_TIPS + 2 <= std::numeric_limits<int16>::max()); // verify int16 ok
    if (a == b) {
      return lg2_rooted[a] + lg2_rooted[n - a];
    }
    // Unify a<b and a>b via lo/hi: removes an unpredictable branch.
    const int16 lo = (a < b) ? a : b;
    const int16 hi = (a < b) ? b : a;
    return lg2_rooted[hi] + lg2_rooted[n - lo] - lg2_rooted[hi - lo + 1];
  }
  
  [[nodiscard]] inline double one_overlap_notb(const int16 a, const int16 n_minus_b, const int16 n) noexcept {
    assert(SL_MAX_TIPS + 2 <= std::numeric_limits<int16>::max()); // verify int16 ok
    const int16 b = n - n_minus_b;
    if (a == b) {
      return lg2_rooted[b] + lg2_rooted[n_minus_b];
    } else if (a < b) {
      return lg2_rooted[b] + lg2_rooted[n - a] - lg2_rooted[b - a + 1];
    } else {
      return lg2_rooted[a] + lg2_rooted[n_minus_b] - lg2_rooted[a - b + 1];
    }
  }


[[nodiscard]] inline double spi_overlap(const splitbit* a_state, const splitbit* b_state,
                     const int16 n_tips, const int16 in_a,
                     const int16 in_b, const int16 n_bins) noexcept {
    
    assert(SL_MAX_BINS <= INT16_MAX);
    
    const splitbit* a_ptr = a_state;
    const splitbit* b_ptr = b_state;
    const splitbit* end_ptr = a_state + n_bins;
    
    bool a_and_b = false;
    
    while(a_ptr != end_ptr) {
      if (*a_ptr & *b_ptr) {
        a_and_b = true;
        break;
      }
      ++a_ptr;
      ++b_ptr;
    }
    
    if (!a_and_b) return one_overlap_notb(in_a, in_b, n_tips);
    
    
    a_ptr = a_state;
    b_ptr = b_state;
    
    bool b_only = false;
    
    while (a_ptr != end_ptr) {
      if (~(*a_ptr) & *b_ptr) {
        b_only = true;
        break;
      }
      ++a_ptr;
      ++b_ptr;
    }
    
    if (!b_only) return one_overlap(in_a, in_b, n_tips);
    
    
    a_ptr = a_state;
    b_ptr = b_state;
    bool a_only = false;
    
    while (a_ptr != end_ptr) {
      if (*a_ptr & ~(*b_ptr)) {
        a_only = true;
        break;
      }
      ++a_ptr;
      ++b_ptr;
    }
    
    if (!a_only) return one_overlap(in_a, in_b, n_tips);
    
    
    const int16 loose_end_tips = n_tips % SL_BIN_SIZE;
    const splitbit tidy_ends = ~(ALL_ONES << loose_end_tips);
    bool neither = false;
    
    for (int16 bin = 0; bin != n_bins; bin++) {
      
      splitbit test = ~(a_state[bin] | b_state[bin]);
      
      if (bin == n_bins - 1 && loose_end_tips) {
        test &= tidy_ends;
      }
      
      if (test) {
        neither = true;
        break;
      }
    }
    
    if (!neither) return one_overlap_notb(in_a, in_b, n_tips);
    
    return 0;
  }
}

extern Rcpp::List cpp_robinson_foulds_distance(Rcpp::RawMatrix x,
                                               Rcpp::RawMatrix y,
                                               Rcpp::IntegerVector nTip);

#endif
