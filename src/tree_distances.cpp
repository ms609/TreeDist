#include <TreeTools/SplitList.h>
#include <TreeDist/popcount.h>
#include <TreeTools/assert.h>
#include <algorithm>
#include <cmath>
#include <memory> /* for unique_ptr, make_unique */
#include <Rcpp/Lightest>
#include "tree_distances.h"
#include "col_prices.h"
#include "split_counts.h"

using Rcpp::_;
using Rcpp::IntegerVector;
using Rcpp::List;
using Rcpp::LogicalVector;
using Rcpp::Named;
using Rcpp::NumericVector;
using Rcpp::RawMatrix;
using TreeTools::SplitList;

namespace TreeDist {
  
  template <typename T>
  inline void resize_uninitialized(std::vector<T>& v, std::size_t n) {
    static_assert(std::is_trivial<T>::value, "Requires trivial type");
    if (n > v.size()) {
      v.reserve(n);
      v.insert(v.end(), n - v.size(), T{});
    } else {
      v.resize(n);
    }
  }

  // Hard C++ ceiling: n_tips^2 must fit int32 and all counters must fit split_int.
  // 32767 = INT16_MAX is the chosen practical cap (n^2 <= 1.07e9 < INT32_MAX).
  // R-level .CheckMaxTips() fires before this; this is a safety backstop.
  static constexpr int32 TREEDIST_MAX_TIPS = 32767;

  void check_ntip(const int32 n) {
    if (n > TREEDIST_MAX_TIPS) {
      Rcpp::stop("Trees with %d tips are not yet supported (maximum %d).",
                 static_cast<int>(n), static_cast<int>(TREEDIST_MAX_TIPS));
    }
  }


}

// [[Rcpp::export]]
int cpp_sl_max_tips() {
  return static_cast<int>(SL_MAX_TIPS);
}

using TreeDist::resize_uninitialized;

inline List robinson_foulds_distance(const RawMatrix &x, const RawMatrix &y,
                                     const int32 n_tips) {

  const SplitList a(x), b(y);
  const int32 last_bin = a.n_bins - 1;
  const int32 unset_tips = (n_tips % SL_BIN_SIZE) ?
    SL_BIN_SIZE - n_tips % SL_BIN_SIZE : 0;
  const splitbit unset_mask = ALL_ONES >> unset_tips;
  cost score = 0;
  
  grf_match matching(a.n_splits, NA_INTEGER);
  
  const int32 n_bins = a.n_bins;
  std::vector<splitbit> b_complement(b.n_splits * n_bins);
  for (int32 i = b.n_splits; i--; ) {
    splitbit* bc_i = &b_complement[i * n_bins];
    for (int32 bin = last_bin; bin--; ) {
      bc_i[bin] = ~b.state[i][bin];
    }
    bc_i[last_bin] = b.state[i][last_bin] ^ unset_mask;
  }
  
  for (int32 ai = a.n_splits; ai--; ) {
    if ((ai & 1023) == 0) Rcpp::checkUserInterrupt();
    for (int32 bi = b.n_splits; bi--; ) {
      
      bool all_match = true;
      bool all_complement = true;
      const splitbit* bc_bi = &b_complement[bi * n_bins];
      
      for (int32 bin = 0; bin < n_bins; ++bin) {
        if ((a.state[ai][bin] != b.state[bi][bin])) {
          all_match = false;
          break;
        }
      }
      if (!all_match) {
        for (int32 bin = 0; bin < n_bins; ++bin) {
          if (a.state[ai][bin] != bc_bi[bin]) {
            all_complement = false;
            break;
          }
        }
      }
      if (all_match || all_complement) {
        ++score;
        matching[ai] = bi + 1;
        break; // Only one match possible per split
      }
    }
  }
  score = cost(a.n_splits + b.n_splits) - score - score;
  
  return List::create(Named("score") = score, _["matching"] = matching);
  
}
template <typename Pop>
inline List robinson_foulds_info_impl(const RawMatrix &x, const RawMatrix &y,
                                 const int32 n_tips) {
  const SplitList a(x), b(y);
  
  const split_int last_bin = a.n_bins - 1;
  const split_int unset_tips = (n_tips % SL_BIN_SIZE) ?
    SL_BIN_SIZE - n_tips % SL_BIN_SIZE : 0;

  const splitbit unset_mask = ALL_ONES >> unset_tips;
  const double lg2_unrooted_n = lg2_unrooted_lookup(n_tips);
  double score = 0;

  grf_match matching(a.n_splits, NA_INTEGER);

  const split_int n_bins = a.n_bins;
  std::vector<splitbit> b_complement(b.n_splits * n_bins);
  for (split_int i = 0; i < b.n_splits; i++) {
    splitbit* bc_i = &b_complement[i * n_bins];
    for (split_int bin = 0; bin < last_bin; ++bin) {
      bc_i[bin] = ~b.state[i][bin];
    }
    bc_i[last_bin] = b.state[i][last_bin] ^ unset_mask;
  }

  for (split_int ai = 0; ai < a.n_splits; ++ai) {
    if ((ai & 1023) == 0) Rcpp::checkUserInterrupt();
    for (split_int bi = 0; bi < b.n_splits; ++bi) {

      bool all_match = true, all_complement = true;
      const splitbit* bc_bi = &b_complement[bi * n_bins];

      for (split_int bin = 0; bin < n_bins; ++bin) {
        if ((a.state[ai][bin] != b.state[bi][bin])) {
          all_match = false;
          break;
        }
      }
      if (!all_match) {
        for (split_int bin = 0; bin < n_bins; ++bin) {
          if ((a.state[ai][bin] != bc_bi[bin])) {
            all_complement = false;
            break;
          }
        }
      }
      if (all_match || all_complement) {
        split_int leaves_in_split = 0;
        for (split_int bin = 0; bin < a.n_bins; ++bin) {
          leaves_in_split += Pop::count(a.state[ai][bin]);
        }

        score += lg2_unrooted_n - lg2_rooted_lookup(leaves_in_split) -
          lg2_rooted_lookup(n_tips - leaves_in_split);

        matching[ai] = bi + 1;
        break; /* Only one match possible per split */
      }
    }
  }
  
  return List::create(Named("score") = score, _["matching"] = matching);
  
}

inline List robinson_foulds_info(const RawMatrix &x, const RawMatrix &y,
                                 const int32 n_tips) {
  return TreeDist::dispatch_popcount([&](auto pop) {
    return robinson_foulds_info_impl<decltype(pop)>(x, y, n_tips);
  });
}



  
template <typename Pop>
inline List matching_split_distance_impl(const RawMatrix &x, const RawMatrix &y,
                                    const int32 n_tips) {
  const SplitList a(x), b(y);
  const split_int most_splits = std::max(a.n_splits, b.n_splits);
  const split_int split_diff = most_splits - std::min(a.n_splits, b.n_splits);
  const split_int half_tips = n_tips / 2;
  if (most_splits == 0) {
    return List::create(Named("score") = 0);
  }
  const cost max_score = BIG / most_splits;
  cost_matrix score(most_splits);

  for (split_int ai = 0; ai < a.n_splits; ++ai) {
    if ((ai & 1023) == 0) Rcpp::checkUserInterrupt();
    for (split_int bi = 0; bi < b.n_splits; ++bi) {
      splitbit total = 0;
      for (split_int bin = 0; bin < a.n_bins; ++bin) {
        total += Pop::count(a.state[ai][bin] ^ b.state[bi][bin]);
      }
      score(ai, bi) = total;
    }
  }

  for (split_int ai = 0; ai < a.n_splits; ++ai) {
    for (split_int bi = 0; bi < b.n_splits; ++bi) {
      if (score(ai, bi) > half_tips) {
        score(ai, bi) = n_tips - score(ai, bi);
      }
    }
    score.padRowAfterCol(ai, b.n_splits, max_score);
  }

  score.padAfterRow(a.n_splits, max_score);

  std::vector<lap_col> rowsol;
  std::vector<lap_row> colsol;

  resize_uninitialized(rowsol, most_splits);
  resize_uninitialized(colsol, most_splits);

  const double final_score = lap(most_splits, score, rowsol, colsol) -
    (max_score * split_diff);

  std::vector<int> final_matching;
  final_matching.reserve(a.n_splits);

  for (split_int i = 0; i < a.n_splits; ++i) {
    const int match = (rowsol[i] < b.n_splits)
    ? static_cast<int>(rowsol[i]) + 1
    : NA_INTEGER;
    final_matching.push_back(match);
  }
  
  
  return List::create(Named("score") = final_score,
                      _["matching"] = final_matching);
}

inline List matching_split_distance(const RawMatrix &x, const RawMatrix &y,
                                    const int32 n_tips) {
  return TreeDist::dispatch_popcount([&](auto pop) {
    return matching_split_distance_impl<decltype(pop)>(x, y, n_tips);
  });
}

template <typename Pop>
inline List jaccard_similarity_impl(const RawMatrix &x, const RawMatrix &y,
                                   const int32 n_tips, const NumericVector &k,
                                   const LogicalVector &allowConflict) {
  const SplitList a(x), b(y);
  const split_int most_splits = std::max(a.n_splits, b.n_splits);

  constexpr cost max_score = BIG;
  constexpr double max_scoreL = max_score;
  const double exponent = k[0];

  bool allow_conflict = allowConflict[0];

  cost_matrix score(most_splits);

  for (split_int ai = 0; ai < a.n_splits; ++ai) {
    if ((ai & 1023) == 0) Rcpp::checkUserInterrupt();

    const split_int na = a.in_split[ai];
    const split_int nA = n_tips - na;

    for (split_int bi = 0; bi < b.n_splits; ++bi) {

      // x divides tips into a|A; y divides tips into b|B
      split_int a_and_b = 0;
      for (split_int bin = 0; bin < a.n_bins; ++bin) {
        a_and_b += Pop::count(a.state[ai][bin] & b.state[bi][bin]);
      }

      const split_int
        nb = b.in_split[bi],
        nB = n_tips - nb,
        a_and_B = na - a_and_b,
        A_and_b = nb - a_and_b,
        A_and_B = nB - a_and_B
      ;

      if (!allow_conflict && !(
            a_and_b == na ||
            a_and_B == na ||
            A_and_b == nA ||
            A_and_B == nA)
          ) {

        score(ai, bi) = max_score; /* Prohibited */

      } else {
        const split_int
          A_or_b  = n_tips - a_and_B,
          a_or_B = n_tips - A_and_b,
          a_or_b = n_tips - A_and_B,
          A_or_B = n_tips - a_and_b
        ;
        const double 
          ars_ab = static_cast<double>(a_and_b) / static_cast<double>(a_or_b),
          ars_Ab = static_cast<double>(A_and_b) / static_cast<double>(A_or_b),
          ars_aB = static_cast<double>(a_and_B) / static_cast<double>(a_or_B),
          ars_AB = static_cast<double>(A_and_B) / static_cast<double>(A_or_B),
          
          min_ars_both = (ars_ab < ars_AB) ? ars_ab : ars_AB,
          min_ars_either = (ars_aB < ars_Ab) ? ars_aB : ars_Ab
        ;
        
        /* LAP will look to minimize an integer. max(ars) is between 0 and 1. */
        if (exponent == 1) {
          /* Nye et al. similarity metric */
          score(ai, bi) = cost(max_scoreL - (max_scoreL * 
            ((min_ars_both > min_ars_either) ? 
            min_ars_both : min_ars_either)));
        } else if (exponent == R_PosInf) {
          score(ai, bi) = cost((min_ars_both == 1 || min_ars_either == 1) ?
                                 0 : max_scoreL);
        } else {
          score(ai, bi) = cost(max_scoreL - (max_scoreL * 
            std::pow((min_ars_both > min_ars_either) ? 
            min_ars_both : min_ars_either, exponent)));
        }
      }
    }
    score.padRowAfterCol(ai, b.n_splits, max_score);
  }
  score.padAfterRow(a.n_splits, max_score);
  
  std::vector<lap_col> rowsol;
  std::vector<lap_row> colsol;
  
  resize_uninitialized(rowsol, most_splits);
  resize_uninitialized(colsol, most_splits);
  
  const double final_score = static_cast<double>((max_score * most_splits) - 
      lap(most_splits, score, rowsol, colsol))
    / max_score;
  
  std::vector<int> final_matching;
  final_matching.reserve(a.n_splits);
  
  for (split_int i = 0; i < a.n_splits; ++i) {
    const int match = (rowsol[i] < b.n_splits)
    ? static_cast<int>(rowsol[i]) + 1
    : NA_INTEGER;
    final_matching.push_back(match);
  }


  return List::create(Named("score") = final_score,
                      _["matching"] = final_matching);

}

inline List jaccard_similarity(const RawMatrix &x, const RawMatrix &y,
                                   const int32 n_tips, const NumericVector &k,
                                   const LogicalVector &allowConflict) {
  return TreeDist::dispatch_popcount([&](auto pop) {
    return jaccard_similarity_impl<decltype(pop)>(x, y, n_tips, k, allowConflict);
  });
}

template <typename Pop>
List msi_distance_impl(const RawMatrix &x, const RawMatrix &y, const int32 n_tips) {
  const SplitList a(x), b(y);
  const split_int most_splits = std::max(a.n_splits, b.n_splits);
  constexpr cost max_score = BIG;
  const double max_possible = lg2_unrooted_lookup(n_tips) -
    lg2_rooted_lookup((n_tips + 1) / 2) - lg2_rooted_lookup(n_tips / 2);
  const double score_over_possible = static_cast<double>(max_score) / max_possible;
  const double possible_over_score = max_possible / max_score;

  cost_matrix score(most_splits);

  // Heap-allocated: n_bins can exceed SL_MAX_BINS for large trees.
  std::vector<splitbit> different(a.n_bins);

  for (split_int ai = 0; ai < a.n_splits; ++ai) {
    if ((ai & 1023) == 0) Rcpp::checkUserInterrupt();
    for (split_int bi = 0; bi < b.n_splits; ++bi) {
      split_int
        n_different = 0,
        n_a_only = 0,
        n_a_and_b = 0
      ;
      for (split_int bin = 0; bin < a.n_bins; ++bin) {
        different[bin] = a.state[ai][bin] ^ b.state[bi][bin];
        n_different += Pop::count(different[bin]);
        n_a_only += Pop::count(a.state[ai][bin] & different[bin]);
        n_a_and_b += Pop::count(a.state[ai][bin] & ~different[bin]);
      }
      const split_int n_same = n_tips - n_different;

      score(ai, bi) = cost(max_score -
        (score_over_possible *
          TreeDist::mmsi_score(n_same, n_a_and_b, n_different, n_a_only)));
    }
    score.padRowAfterCol(ai, b.n_splits, max_score);
  }
  score.padAfterRow(a.n_splits, max_score);

  std::vector<lap_col> rowsol;
  std::vector<lap_row> colsol;

  resize_uninitialized(rowsol, most_splits);
  resize_uninitialized(colsol, most_splits);

  const double final_score = static_cast<double>(
    (max_score * most_splits) - lap(most_splits, score, rowsol, colsol)) *
      possible_over_score;

  std::vector<int> final_matching;
  final_matching.reserve(a.n_splits);

  for (split_int i = 0; i < a.n_splits; ++i) {
    const int match = (rowsol[i] < b.n_splits)
    ? static_cast<int>(rowsol[i]) + 1
    : NA_INTEGER;
    final_matching.push_back(match);
  }

  return List::create(Named("score") = final_score,
                      _["matching"] = final_matching);

}

List msi_distance(const RawMatrix &x, const RawMatrix &y, const int32 n_tips) {
  return TreeDist::dispatch_popcount([&](auto pop) {
    return msi_distance_impl<decltype(pop)>(x, y, n_tips);
  });
}

// Warm-started LAP; see lap_warm.cpp.
extern cost lap_warm(lap_row dim, CostMatrix &input_cost,
                     std::vector<lap_col> &rowsol,
                     std::vector<lap_row> &colsol, LapScratch &scratch,
                     const cost *v_init);

// With Bounded = true, gives up once the score is known to fall below
// min_score.  It then returns an upper bound on the score, below min_score,
// and an empty matching.  Scores at or above min_score are computed exactly as
// with Bounded = false.
//
// Bounded only: `prices` (may be null) can supply a price v_j for each split
// j of `y`, typically the optimal column duals of a closely related problem.
// When both trees have the same number of splits, the matching is a perfect
// assignment, and for any prices
//   sum_j v_j + sum_i min_j (c_ij - v_j)  <=  minimum assignment cost,
// because u_i = min_j (c_ij - v_j) makes (u, v) feasible for the dual.  This
// gives an upper bound on the score once the cost matrix is filled.  If that
// bound does not reject the candidate, the prices warm-start the LAP.
// Prices affect only speed, never the result.
// Reusable cost matrices for mutual_clustering<true>().
static cost_matrix &BoundedPool(const int i, const std::size_t dim) {
  static thread_local cost_matrix pool[2];
  pool[i].resize(dim);
  return pool[i];
}

// k log2 k for 0 <= k <= n_tips, for the bounded mutual clustering fill.
static const double *KLg2K(const int32 n_tips) {
  static thread_local std::vector<double> table;
  if (table.size() <= static_cast<std::size_t>(n_tips)) {
    const int32 old_size = static_cast<int32>(table.size());
    table.resize(n_tips + 1);
    for (int32 k = old_size; k <= n_tips; ++k) {
      table[k] = k * lg2_lookup(k);
    }
  }
  return table.data();
}

template <bool Bounded, typename Pop>
List mutual_clustering_impl(const RawMatrix &x, const RawMatrix &y,
                       const int32 n_tips, const double min_score,
                       const ColPrices *prices) {
  const SplitList a(x);
  const SplitList b(y);
  const bool a_has_more_splits = (a.n_splits > b.n_splits);
  const split_int most_splits = a_has_more_splits ? a.n_splits : b.n_splits;
  const split_int a_extra_splits = a_has_more_splits ? most_splits - b.n_splits : 0;
  const split_int b_extra_splits = a_has_more_splits ? 0 : most_splits - a.n_splits;
  const double n_tips_reciprocal = 1.0 / n_tips;

  if (most_splits == 0 || n_tips == 0) {
    return List::create(Named("score") = 0,
                        _["matching"] = IntegerVector(0));
  }

  constexpr cost max_score = BIG;
  constexpr double over_max_score = 1.0 / static_cast<double>(max_score);
  const double max_over_tips = static_cast<double>(max_score) * n_tips_reciprocal;
  const double lg2_n = lg2_lookup(n_tips);

  // Bounded only: upper bound on the score.  A split's mutual clustering
  // information with any other split cannot exceed its own entropy.
  [[maybe_unused]] double upper_bound = 0;
  if constexpr (Bounded) {
    double entropy_a = 0;
    for (split_int ai = 0; ai < a.n_splits; ++ai) {
      entropy_a += TreeDist::ic_matching(a.in_split[ai],
                                         n_tips - a.in_split[ai], n_tips);
    }
    double entropy_b = 0;
    for (split_int bi = 0; bi < b.n_splits; ++bi) {
      entropy_b += TreeDist::ic_matching(b.in_split[bi],
                                         n_tips - b.in_split[bi], n_tips);
    }
    upper_bound = std::min(entropy_a, entropy_b) * n_tips_reciprocal;
    if (upper_bound < min_score) {
      return List::create(Named("score") = upper_bound,
                          _["matching"] = IntegerVector(0));
    }
    upper_bound = entropy_a * n_tips_reciprocal;
  }

  // Bounded only: called once per candidate in SPR leaf scoring, so reuse
  // buffers rather than reallocating them on each call.
  cost_matrix local_score(Bounded ? 0 : most_splits);
  cost_matrix &score = Bounded ? BoundedPool(0, most_splits) : local_score;

  double exact_match_score = 0;
  split_int exact_matches = 0;
  // vector zero-initializes [so does make_unique]
  // match will have one added to it so numbering follows R; hence 0 = UNMATCHED
  std::vector<int> a_match(a.n_splits);
  std::unique_ptr<split_int[]> b_match = std::make_unique<split_int[]>(b.n_splits);

  // Bounded only: write each cell's information as
  //   sum_q f(n_q) - f(na) - f(nA) - f(nb) - f(nB) + f(n),
  // where f(k) = k log2 k and q ranges over the four quadrants, hoisting the
  // row and column terms.  This takes four table lookups and additions per
  // cell, against four lookups and multiplications; the result agrees with
  // the unbounded form to within a unit of the integer cost.
  [[maybe_unused]] static thread_local std::vector<double> col_term;
  [[maybe_unused]] const double *f = nullptr;
  if constexpr (Bounded) {
    f = KLg2K(n_tips);
    col_term.resize(b.n_splits);
    for (split_int bi = 0; bi < b.n_splits; ++bi) {
      col_term[bi] = -f[b.in_split[bi]] - f[n_tips - b.in_split[bi]];
    }
  }

  for (split_int ai = 0; ai < a.n_splits; ++ai) {
    if (!Bounded && (ai & 1023) == 0) Rcpp::checkUserInterrupt();
    if (a_match[ai]) continue;

    const split_int na = a.in_split[ai];
    const split_int nA = n_tips - na;
    const auto *a_row = a.state[ai];
    [[maybe_unused]] const double offset_a = lg2_n - lg2_lookup(na);
    [[maybe_unused]] const double offset_A = lg2_n - lg2_lookup(nA);
    [[maybe_unused]] cost row_min = max_score;
    [[maybe_unused]] double row_term = 0;
    if constexpr (Bounded) {
      row_term = f[n_tips] - f[na] - f[nA];
    }

    for (split_int bi = 0; bi < b.n_splits; ++bi) {

      // x divides tips into a|A; y divides tips into b|B
      split_int a_and_b = 0;
      const auto *b_row = b.state[bi];
      for (split_int bin = 0; bin < a.n_bins; ++bin) {
        a_and_b += Pop::count(a_row[bin] & b_row[bin]);
      }

      const split_int nb = b.in_split[bi];
      const split_int nB = n_tips - nb;
      const split_int a_and_B = na - a_and_b;
      const split_int A_and_b = nb - a_and_b;
      const split_int A_and_B = nA - A_and_b;

      if ((!a_and_B && !A_and_b) ||
          (!a_and_b && !A_and_B)) {
        exact_match_score += TreeDist::ic_matching(na, nA, n_tips);
        ++exact_matches;
        a_match[ai] = bi + 1;
        b_match[bi] = ai + 1;
        break;
      } else if (a_and_b == A_and_b
                   && a_and_b == a_and_B
                   && a_and_b == A_and_B) {
        score(ai, bi) = max_score; // Avoid rounding errors
      } else if constexpr (Bounded) {
        const double ic_sum = f[a_and_b] + f[a_and_B] + f[A_and_b] +
          f[A_and_B] + row_term + col_term[bi];
        score(ai, bi) = max_score - static_cast<cost>(ic_sum * max_over_tips);
      } else {
        const double lg2_nb = lg2_lookup(nb);
        const double lg2_nB = lg2_lookup(nB);
        const double ic_sum =
          a_and_b * (lg2_lookup(a_and_b) + offset_a - lg2_nb) +
          a_and_B * (lg2_lookup(a_and_B) + offset_a - lg2_nB) +
          A_and_b * (lg2_lookup(A_and_b) + offset_A - lg2_nb) +
          A_and_B * (lg2_lookup(A_and_B) + offset_A - lg2_nB);

        // Division by n_tips converts n(A&B) to P(A&B) for each ic_element
        score(ai, bi) = max_score - static_cast<cost>(ic_sum * max_over_tips);
      }
      if constexpr (Bounded) {
        if (score(ai, bi) < row_min) row_min = score(ai, bi);
      }
    }

    if constexpr (Bounded) {
      if (!a_match[ai]) {
        // Replace this split's entropy with its best available match.
        upper_bound += static_cast<double>(max_score - row_min) *
          over_max_score - TreeDist::ic_matching(na, nA, n_tips) *
          n_tips_reciprocal;
        if (upper_bound < min_score) {
          return List::create(Named("score") = upper_bound,
                              _["matching"] = IntegerVector(0));
        }
      }
    }

    if (b.n_splits < most_splits) {
      score.padRowAfterCol(ai, b.n_splits, max_score);
    }
  }

  if (exact_matches == b.n_splits || exact_matches == a.n_splits) {
    return List::create(
      Named("score") = exact_match_score * n_tips_reciprocal,
      _["matching"] = a_match);
  }

  const split_int lap_dim = most_splits - exact_matches;
  ASSERT(lap_dim > 0);

  // Bounded: dual bound from prices, on the problem that remains once exact
  // matches are removed.  Each row minimum ranges over every column, so it
  // cannot exceed the minimum over the remaining columns.  Prices lie in
  // [0, BIG] and costs in [0, BIG], so |bound| <= 2 * lap_dim * BIG, which
  // fits in a cost when most_splits <= SL_MAX_SPLITS / 2.
  [[maybe_unused]] const cost *col_price = nullptr;
  [[maybe_unused]] static thread_local std::vector<cost> v_init;
  if constexpr (Bounded) {
    if (prices && a.n_splits == b.n_splits &&
        most_splits <= SL_MAX_SPLITS / 2) {
      col_price = prices->get(prices->context);
    }
    if (col_price) {
      cost lb = 0;
      v_init.clear();
      for (split_int bi = 0; bi < b.n_splits; ++bi) {
        if (!b_match[bi]) {
          lb += col_price[bi];
          v_init.push_back(col_price[bi]);
        }
      }
      for (split_int ai = 0; ai < a.n_splits; ++ai) {
        if (a_match[ai]) continue;
        const cost *row = score.row(ai);
        cost row_min = row[0] - col_price[0];
        for (split_int bi = 1; bi < b.n_splits; ++bi) {
          const cost reduced = row[bi] - col_price[bi];
          if (reduced < row_min) row_min = reduced;
        }
        lb += row_min;
      }
      const cost lap_ub = max_score * lap_dim - lb;
      const double bound = exact_matches ?
        static_cast<double>(lap_ub) * over_max_score +
          exact_match_score / n_tips :
        static_cast<double>(lap_ub) / max_score;
      if (bound < min_score) {
        return List::create(Named("score") = bound,
                            _["matching"] = IntegerVector(0));
      }
    }
  }

  std::vector<lap_col> rowsol;
  std::vector<lap_row> colsol;
  resize_uninitialized(rowsol, lap_dim);
  resize_uninitialized(colsol, lap_dim);
  [[maybe_unused]] static thread_local LapScratch warm_scratch;
  const auto Solve = [&](cost_matrix &m) -> cost {
    if constexpr (Bounded) {
      if (col_price) {
        return lap_warm(lap_dim, m, rowsol, colsol, warm_scratch,
                        v_init.data());
      }
    }
    return lap(lap_dim, m, rowsol, colsol, !Bounded);
  };
  cost_matrix local_small(Bounded ? 0 : lap_dim);
  cost_matrix &small_score = Bounded ? BoundedPool(1, lap_dim) : local_small;

  if (exact_matches) {
    split_int a_pos = 0;
    for (split_int ai = 0; ai < a.n_splits; ++ai) {
      if (a_match[ai]) continue;
      split_int b_pos = 0;
      for (split_int bi = 0; bi < b.n_splits; ++bi) {
        if (b_match[bi]) continue;
        small_score(a_pos, b_pos) = score(ai, bi);
        b_pos++;
      }
      small_score.padRowAfterCol(a_pos, lap_dim - a_extra_splits, max_score);
      a_pos++;
    }
    small_score.padAfterRow(lap_dim - b_extra_splits, max_score);

    const double lap_score = static_cast<double>((max_score * lap_dim) -
      Solve(small_score)) * over_max_score;
    const double final_score = lap_score + (exact_match_score / n_tips);

    std::unique_ptr<split_int[]> lap_decode = std::make_unique<split_int[]>(lap_dim);
    split_int fuzzy_match = 0;
    for (split_int bi = 0; bi < b.n_splits; ++bi) {
      if (!b_match[bi]) {
        assert(fuzzy_match < lap_dim);
        lap_decode[fuzzy_match++] = bi + 1;
      }
    }

    fuzzy_match = 0;
    std::vector<int> final_matching;
    TreeDist::resize_uninitialized(final_matching, a.n_splits);
    for (split_int i = 0; i < a.n_splits; ++i) {
      if (a_match[i]) {
        final_matching[i] = a_match[i];
      } else {
        assert(fuzzy_match < lap_dim);
        const split_int row_idx = fuzzy_match++;
        const split_int col_idx = rowsol[row_idx];
        final_matching[i] = (col_idx >= lap_dim - a_extra_splits) ? NA_INTEGER :
          lap_decode[col_idx];
      }
    }

    return List::create(Named("score") = final_score,
                        _["matching"] = final_matching);
  } else {

    for (split_int ai = a.n_splits; ai < most_splits; ++ai) {
      for (split_int bi = 0; bi < most_splits; ++bi) {
        score(ai, bi) = max_score;
      }
    }

    const double final_score = static_cast<double>(
      (max_score * lap_dim) - Solve(score)
    ) / max_score;

    std::vector<int> final_matching;
    final_matching.reserve(a.n_splits);
    for (split_int i = 0; i < a.n_splits; ++i) {
      const int match = (rowsol[i] < b.n_splits)
      ? static_cast<int>(rowsol[i]) + 1
      : NA_INTEGER;
      final_matching.push_back(match);
    }

    return List::create(Named("score") = final_score,
                        _["matching"] = final_matching);
  }
}

template <bool Bounded>
List mutual_clustering(const RawMatrix &x, const RawMatrix &y,
                       const int32 n_tips, const double min_score = 0,
                       const ColPrices *prices = nullptr) {
  return TreeDist::dispatch_popcount([&](auto pop) {
    return mutual_clustering_impl<Bounded, decltype(pop)>(x, y, n_tips,
                                                          min_score, prices);
  });
}

template <bool Fast, typename Pop>
static inline void shared_phylo_fill_impl(const SplitList& a, const SplitList& b,
                                     const int32 n_tips,
                                     const double best_overlap,
                                     const double score_over_possible,
                                     const cost max_score,
                                     cost_matrix& score) {
  for (split_int ai = a.n_splits; ai--; ) {
    if ((ai & 1023) == 0) Rcpp::checkUserInterrupt();
    for (split_int bi = b.n_splits; bi--; ) {
      const double spi_over = TreeDist::spi_overlap<Fast, Pop>(
        a.state[ai], b.state[bi], n_tips, a.in_split[ai], b.in_split[bi],
                                                                    a.n_bins);

      score(ai, bi) = spi_over < 0.0 ? max_score :
        cost((spi_over - best_overlap) * score_over_possible);
    }
    score.padRowAfterCol(ai, b.n_splits, max_score);
  }
}

// Templated LAP-fill so the per-cell spi_overlap resolves to direct
// lg2_rooted[] access on the hot path (n_tips <= SL_MAX_TIPS).
template<bool Fast>
static inline void shared_phylo_fill(const SplitList& a, const SplitList& b,
                                     const int32 n_tips,
                                     const double best_overlap,
                                     const double score_over_possible,
                                     const cost max_score,
                                     cost_matrix& score) {
  return TreeDist::dispatch_popcount([&](auto pop) {
    return shared_phylo_fill_impl<Fast, decltype(pop)>(a, b, n_tips, best_overlap, score_over_possible, max_score, score);
  });
}

inline List shared_phylo (const RawMatrix &x, const RawMatrix &y,
                          const int32 n_tips,
                          const bool force_slow = false) {
  const SplitList a(x), b(y);
  const split_int most_splits = std::max(a.n_splits, b.n_splits);
  const split_int overlap_a = (n_tips + 1) / 2;

  constexpr cost max_score = BIG;
  // force_slow: test-only hook so the bounds-checked lookup path can be
  // exercised on small trees for code coverage of the template <false>
  // instantiations and the lg2_*_slow fallbacks.
  const bool fast = !force_slow && n_tips <= static_cast<int32>(SL_MAX_TIPS + 1);
  double lg2_unrooted_n, best_overlap;
  if (fast) {
    lg2_unrooted_n = lg2_unrooted[n_tips];
    best_overlap = TreeDist::one_overlap<true>(overlap_a, n_tips / 2, n_tips);
  } else {
    lg2_unrooted_n = lg2_unrooted_lookup(n_tips);
    best_overlap = TreeDist::one_overlap<false>(overlap_a, n_tips / 2, n_tips);
  }
  const double max_possible = lg2_unrooted_n - best_overlap;
  const double score_over_possible = max_score / max_possible;
  const double possible_over_score = max_possible / max_score;

  // a and b are "clades" separating an "ingroup" [1] from an "outgroup" [0].
  // In/out direction [i.e. 1/0 bit] is arbitrary.
  cost_matrix score(most_splits);

  if (fast) {
    shared_phylo_fill<true>(a, b, n_tips, best_overlap,
                            score_over_possible, max_score, score);
  } else {
    shared_phylo_fill<false>(a, b, n_tips, best_overlap,
                             score_over_possible, max_score, score);
  }
  score.padAfterRow(a.n_splits, max_score);

  std::vector<lap_col> rowsol;
  std::vector<lap_row> colsol;

  resize_uninitialized(rowsol, most_splits);
  resize_uninitialized(colsol, most_splits);

  const double final_score = static_cast<double>(
      (max_score * most_splits) - lap(most_splits, score, rowsol, colsol)) *
        possible_over_score;

  std::vector<int> final_matching;
  final_matching.reserve(a.n_splits);

  for (split_int i = 0; i < a.n_splits; ++i) {
    const int match = (rowsol[i] < b.n_splits)
    ? static_cast<int>(rowsol[i]) + 1
    : NA_INTEGER;
    final_matching.push_back(match);
  }

  return List::create(Named("score") = final_score,
                      _["matching"] = final_matching);
}

// [[Rcpp::export]]
List cpp_robinson_foulds_distance(const RawMatrix &x, const RawMatrix &y,
                                  const IntegerVector &nTip) {
  ASSERT(x.cols() == y.cols() && "Input splits must address same number of tips.");
  const int32 n_tip = static_cast<int32>(nTip[0]);
  TreeDist::check_ntip(n_tip);
  return robinson_foulds_distance(x, y, n_tip);
}

// [[Rcpp::export]]
List cpp_robinson_foulds_info(const RawMatrix &x, const RawMatrix &y,
                              const IntegerVector &nTip) {
  ASSERT(x.cols() == y.cols() && "Input splits must address same number of tips.");
  const int32 n_tip = static_cast<int32>(nTip[0]);
  TreeDist::check_ntip(n_tip);
  return robinson_foulds_info(x, y, n_tip);
}

// [[Rcpp::export]]
List cpp_matching_split_distance(const RawMatrix &x, const RawMatrix &y,
                                 const IntegerVector &nTip) {
  ASSERT(x.cols() == y.cols() && "Input splits must address same number of tips.");
  const int32 n_tip = static_cast<int32>(nTip[0]);
  TreeDist::check_ntip(n_tip);
  return matching_split_distance(x, y, n_tip);
}

// [[Rcpp::export]]
List cpp_jaccard_similarity(const RawMatrix &x, const RawMatrix &y,
                            const IntegerVector &nTip, const NumericVector &k,
                            const LogicalVector &allowConflict) {
  ASSERT(x.cols() == y.cols() && "Input splits must address same number of tips.");
  const int32 n_tip = static_cast<int32>(nTip[0]);
  TreeDist::check_ntip(n_tip);
  return jaccard_similarity(x, y, n_tip, k, allowConflict);
}

// [[Rcpp::export]]
List cpp_msi_distance(const RawMatrix &x, const RawMatrix &y,
                      const IntegerVector &nTip) {
  ASSERT(x.cols() == y.cols() && "Input splits must address same number of tips.");
  const int32 n_tip = static_cast<int32>(nTip[0]);
  TreeDist::check_ntip(n_tip);
  return msi_distance(x, y, n_tip);
}

// [[Rcpp::export]]
List cpp_mutual_clustering(const RawMatrix &x, const RawMatrix &y,
                           const IntegerVector &nTip) {
  ASSERT(x.cols() == y.cols() && "Input splits must address same number of tips.");
  const int32 n_tip = static_cast<int32>(nTip[0]);
  TreeDist::check_ntip(n_tip);
  return mutual_clustering<false>(x, y, n_tip);
}

// Mutual clustering information, or an upper bound on it if it is below
// min_score.  Used to discard candidates early in SPRDist(method = "Rogue").
double mutual_clustering_bounded_score(const RawMatrix &x, const RawMatrix &y,
                                       const int32 n_tip,
                                       const double min_score,
                                       const ColPrices *prices) {
  return Rcpp::as<double>(
    mutual_clustering<true>(x, y, n_tip, min_score, prices)["score"]);
}

// As mutual_clustering_bounded_score(), for trees obtained by removing one
// or two leaves from a parent pair whose splits overlap by `sc.overlap`, when
// the resulting trees share no split.  See split_counts.h.
template <bool Pair>
double bounded_counts_impl(const SplitCounts &sc, const double min_score,
                           const ColPrices *prices) {
  const split_int n_tips = sc.n_tip;
  const split_int dim = sc.n_split;
  constexpr cost max_score = BIG;
  constexpr double over_max_score = 1.0 / static_cast<double>(max_score);
  const double n_tips_reciprocal = 1.0 / n_tips;
  const double max_over_tips = static_cast<double>(max_score) *
    n_tips_reciprocal;

  // Row bound, as mutual_clustering<true>().  The entropy bound that precedes
  // it there is omitted: it rejects a candidate only if the trees' entropies
  // differ by more than the best score so far, which no candidate met in
  // 600,000 tested.
  double upper_bound = 0;
  for (split_int i = 0; i < dim; ++i) {
    upper_bound += TreeDist::ic_matching(sc.row_in[i], n_tips - sc.row_in[i],
                                         n_tips);
  }
  upper_bound *= n_tips_reciprocal;

  const double *f = KLg2K(n_tips);
  static thread_local std::vector<double> col_term;
  col_term.resize(dim);
  for (split_int j = 0; j < dim; ++j) {
    col_term[j] = -f[sc.col_in[j]] - f[n_tips - sc.col_in[j]];
  }
  cost_matrix &score = BoundedPool(0, dim);
  for (split_int i = 0; i < dim; ++i) {
    const split_int na = sc.row_in[i];
    const split_int nA = n_tips - na;
    const int16_t *overlap = sc.overlap + std::size_t(sc.row[i]) * sc.stride;
    const int16_t row_has = sc.row_has[i];
    const double row_term = f[n_tips] - f[na] - f[nA];
    cost *out = score.row(i);
    cost row_min = max_score;
    for (split_int j = 0; j < dim; ++j) {
      const int16_t both = row_has & sc.col_has[j];
      const split_int a_and_b = overlap[sc.col[j]] -
        (Pair ? (both & 1) + (both >> 1) : both);
      const split_int nb = sc.col_in[j];
      const split_int a_and_B = na - a_and_b;
      const split_int A_and_b = nb - a_and_b;
      const split_int A_and_B = nA - A_and_b;
      ASSERT(a_and_B || A_and_b);
      ASSERT(a_and_b || A_and_B);
      const cost c = (a_and_b == A_and_b && a_and_b == a_and_B &&
                      a_and_b == A_and_B) ? max_score :
        max_score - static_cast<cost>(
            (f[a_and_b] + f[a_and_B] + f[A_and_b] + f[A_and_B] + row_term +
            col_term[j]) * max_over_tips);
      out[j] = c;
      if (c < row_min) row_min = c;
    }
    upper_bound += static_cast<double>(max_score - row_min) * over_max_score -
      TreeDist::ic_matching(na, nA, n_tips) * n_tips_reciprocal;
    if (upper_bound < min_score) {
      return upper_bound;
    }
  }

  const cost *col_price = (prices && dim <= SL_MAX_SPLITS / 2) ?
    prices->get(prices->context) : nullptr;
  if (col_price) {
    cost lb = 0;
    for (split_int j = 0; j < dim; ++j) {
      lb += col_price[j];
    }
    for (split_int i = 0; i < dim; ++i) {
      const cost *row = score.row(i);
      cost row_min = row[0] - col_price[0];
      for (split_int j = 1; j < dim; ++j) {
        const cost reduced = row[j] - col_price[j];
        if (reduced < row_min) row_min = reduced;
      }
      lb += row_min;
    }
    const double bound = static_cast<double>(max_score * dim - lb) *
      over_max_score;
    if (bound < min_score) {
      return bound;
    }
  }

  static thread_local std::vector<lap_col> rowsol;
  static thread_local std::vector<lap_row> colsol;
  static thread_local LapScratch warm_scratch;
  resize_uninitialized(rowsol, dim);
  resize_uninitialized(colsol, dim);
  const cost lap_cost = col_price ?
    lap_warm(dim, score, rowsol, colsol, warm_scratch, col_price) :
    lap(dim, score, rowsol, colsol, false); // #nocov: > 16k leaves
  return static_cast<double>(max_score * dim - lap_cost) * over_max_score;
}

double mutual_clustering_bounded_counts(const SplitCounts &sc,
                                        const double min_score,
                                        const ColPrices *prices) {
  return sc.pair ? bounded_counts_impl<true>(sc, min_score, prices) :
    bounded_counts_impl<false>(sc, min_score, prices);
}

// Optimal column duals of the mutual clustering LAP between the splits `x`
// and `y` (identical splits included), one per split in `y`, in cost units,
// shifted so that the least is zero and capped at BIG: prices for
// mutual_clustering<true>() on related problems.  Empty if there is nothing
// to match.
std::vector<cost> mutual_clustering_col_prices(const RawMatrix &x,
                                               const RawMatrix &y,
                                               const int32 n_tips) {
  const SplitList a(x);
  const SplitList b(y);
  const split_int dim = std::max(a.n_splits, b.n_splits);
  if (dim == 0 || n_tips < 4 || b.n_splits == 0) {
    return {}; // #nocov
  }
  constexpr cost max_score = BIG;
  const double max_over_tips = static_cast<double>(max_score) / n_tips;
  const double lg2_n = lg2_lookup(n_tips);
  cost_matrix score(dim);
  for (split_int ai = 0; ai < a.n_splits; ++ai) {
    const split_int na = a.in_split[ai];
    const split_int nA = n_tips - na;
    const double offset_a = lg2_n - lg2_lookup(na);
    const double offset_A = lg2_n - lg2_lookup(nA);
    for (split_int bi = 0; bi < b.n_splits; ++bi) {
      split_int a_and_b = 0;
      for (split_int bin = 0; bin < a.n_bins; ++bin) {
        a_and_b += TreeDist::count_bits(a.state[ai][bin] & b.state[bi][bin]);
      }
      const split_int nb = b.in_split[bi];
      const split_int nB = n_tips - nb;
      const split_int a_and_B = na - a_and_b;
      const split_int A_and_b = nb - a_and_b;
      const split_int A_and_B = nA - A_and_b;
      if (a_and_b == A_and_b && a_and_b == a_and_B && a_and_b == A_and_B) {
        score(ai, bi) = max_score;
      } else {
        const double lg2_nb = lg2_lookup(nb);
        const double lg2_nB = lg2_lookup(nB);
        const double ic_sum =
          a_and_b * (lg2_lookup(a_and_b) + offset_a - lg2_nb) +
          a_and_B * (lg2_lookup(a_and_B) + offset_a - lg2_nB) +
          A_and_b * (lg2_lookup(A_and_b) + offset_A - lg2_nb) +
          A_and_B * (lg2_lookup(A_and_B) + offset_A - lg2_nB);
        score(ai, bi) = max_score - static_cast<cost>(ic_sum * max_over_tips);
      }
    }
    score.padRowAfterCol(ai, b.n_splits, max_score);
  }
  score.padAfterRow(a.n_splits, max_score);
  std::vector<lap_col> rowsol(dim);
  std::vector<lap_row> colsol(dim);
  LapScratch scratch;
  lap(dim, score, rowsol, colsol, false, scratch);
  std::vector<cost> ret(scratch.v.begin(), scratch.v.begin() + b.n_splits);
  const cost least = *std::min_element(ret.begin(), ret.end());
  for (cost &v : ret) {
    v = std::min(v - least, max_score);
  }
  return ret;
}

// [[Rcpp::export]]
double cpp_mutual_clustering_bounded(const RawMatrix &x, const RawMatrix &y,
                                     const IntegerVector &nTip,
                                     const double minScore) {
  ASSERT(x.cols() == y.cols() && "Input splits must address same number of tips.");
  const int32 n_tip = static_cast<int32>(nTip[0]);
  TreeDist::check_ntip(n_tip);
  return mutual_clustering_bounded_score(x, y, n_tip, minScore, nullptr);
}

// [[Rcpp::export]]
List cpp_shared_phylo(const RawMatrix &x, const RawMatrix &y,
                      const IntegerVector &nTip,
                      const bool force_slow = false) {
  ASSERT(x.cols() == y.cols() && "Input splits must address same number of tips.");
  const int32 n_tip = static_cast<int32>(nTip[0]);
  TreeDist::check_ntip(n_tip);
  return shared_phylo(x, y, n_tip, force_slow);
}

// Test hook: force the portable popcount path so that tests exercise the
// code that runs on CPUs without POPCNT. Returns the previous setting.
// [[Rcpp::export]]
bool cpp_popcount_force_sw(const bool force) {
  const bool previous = TreeDist::popcount_force_sw;
  TreeDist::popcount_force_sw = force;
  return previous;
}
