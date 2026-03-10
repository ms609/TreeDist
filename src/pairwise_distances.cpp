/* pairwise_distances.cpp
 *
 * Batch pairwise distance functions using OpenMP where available.
 * Each exported function takes a list of split matrices (one per tree) and
 * returns the lower-triangle distance vector in combn(n, 2) order, suitable
 * for direct use as the data payload of an R dist object.
 */

#ifdef _OPENMP
#include <omp.h>
#endif

#include <TreeTools/SplitList.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <vector>
#include <Rcpp/Lightest>

#include "tree_distances.h"

using Rcpp::List;
using Rcpp::NumericVector;
using Rcpp::RawMatrix;
using TreeTools::SplitList;
using TreeTools::count_bits;


// Score-only version of mutual_clustering(). Thread-safe: uses only local
// storage and read-only globals (lg2 table, lookup tables).
// Passes allow_interrupt = false to lap() so it is safe to call from an
// OpenMP parallel region.
static double mutual_clustering_score(
    const SplitList& a, const SplitList& b, const int32 n_tips
) {
  if (a.n_splits == 0 || b.n_splits == 0 || n_tips == 0) return 0.0;

  const bool a_has_more = (a.n_splits > b.n_splits);
  const int16 most_splits = a_has_more ? a.n_splits : b.n_splits;
  const int16 a_extra     = a_has_more ? most_splits - b.n_splits : 0;
  const int16 b_extra     = a_has_more ? 0 : most_splits - a.n_splits;
  const double n_tips_rcp = 1.0 / static_cast<double>(n_tips);

  constexpr cost max_score  = BIG;
  constexpr double over_max = 1.0 / static_cast<double>(BIG);
  const double max_over_tips = static_cast<double>(BIG) * n_tips_rcp;

  cost_matrix score(most_splits);
  double exact_score  = 0.0;
  int16  exact_n      = 0;

  std::vector<int>  a_match(a.n_splits, 0);
  auto b_match = std::make_unique<int16[]>(b.n_splits);
  std::fill(b_match.get(), b_match.get() + b.n_splits, int16(0));

  for (int16 ai = 0; ai < a.n_splits; ++ai) {
    if (a_match[ai]) continue;

    const int16  na    = a.in_split[ai];
    const int16  nA    = n_tips - na;
    const auto*  a_row = a.state[ai];

    for (int16 bi = 0; bi < b.n_splits; ++bi) {
      int16       a_and_b = 0;
      const auto* b_row   = b.state[bi];
      for (int16 bin = 0; bin < a.n_bins; ++bin) {
        a_and_b += count_bits(a_row[bin] & b_row[bin]);
      }

      const int16 nb    = b.in_split[bi];
      const int16 nB    = n_tips - nb;
      const int16 a_and_B = na - a_and_b;
      const int16 A_and_b = nb - a_and_b;
      const int16 A_and_B = nA - A_and_b;

      if ((!a_and_B && !A_and_b) || (!a_and_b && !A_and_B)) {
        // Exact match (nested or identical splits)
        exact_score += TreeDist::ic_matching(na, nA, n_tips);
        ++exact_n;
        a_match[ai] = bi + 1;
        b_match[bi] = ai + 1;
        break;
      } else if (a_and_b == A_and_b && a_and_b == a_and_B && a_and_b == A_and_B) {
        score(ai, bi) = max_score; // Avoid rounding errors on orthogonal splits
      } else {
        double ic_sum = 0.0;
        TreeDist::add_ic_element(ic_sum, a_and_b, na, nb, n_tips);
        TreeDist::add_ic_element(ic_sum, a_and_B, na, nB, n_tips);
        TreeDist::add_ic_element(ic_sum, A_and_b, nA, nb, n_tips);
        TreeDist::add_ic_element(ic_sum, A_and_B, nA, nB, n_tips);
        score(ai, bi) = max_score - static_cast<cost>(ic_sum * max_over_tips);
      }
    }

    if (b.n_splits < most_splits) {
      score.padRowAfterCol(ai, b.n_splits, max_score);
    }
  }

  // Early exit when everything matched exactly
  if (exact_n == b.n_splits || exact_n == a.n_splits) {
    return exact_score * n_tips_rcp;
  }

  const int16 lap_n = most_splits - exact_n;
  std::vector<lap_col> rowsol(lap_n);
  std::vector<lap_row> colsol(lap_n);

  if (exact_n) {
    // Build a reduced cost matrix omitting exact-matched rows/cols
    cost_matrix small(lap_n);
    int16 a_pos = 0;
    for (int16 ai = 0; ai < a.n_splits; ++ai) {
      if (a_match[ai]) continue;
      int16 b_pos = 0;
      for (int16 bi = 0; bi < b.n_splits; ++bi) {
        if (b_match[bi]) continue;
        small(a_pos, b_pos) = score(ai, bi);
        ++b_pos;
      }
      small.padRowAfterCol(a_pos, lap_n - a_extra, max_score);
      ++a_pos;
    }
    small.padAfterRow(lap_n - b_extra, max_score);

    const double lap_score =
      static_cast<double>((max_score * lap_n) -
                          lap(lap_n, small, rowsol, colsol, false)) * over_max;
    return lap_score + exact_score * n_tips_rcp;

  } else {
    // No exact matches — pad and solve the full matrix
    for (int16 ai = a.n_splits; ai < most_splits; ++ai) {
      for (int16 bi = 0; bi < most_splits; ++bi) {
        score(ai, bi) = max_score;
      }
    }
    return static_cast<double>(
      (max_score * lap_n) - lap(lap_n, score, rowsol, colsol, false)
    ) / max_score;
  }
}


//' Pairwise mutual clustering information — batch computation
//'
//' Internal function. Computes all pairwise MCI scores for a set of trees,
//' using OpenMP threads when available (falling back to single-threaded
//' execution otherwise). No interrupt checking is performed inside the
//' parallel region; the outer R call remains interruptible between batches.
//'
//' @param splits_list A list of split matrices (class `Splits` or `RawMatrix`),
//'   one per tree, all covering the same tip set.  Typically the object
//'   returned by `as.Splits(trees, tipLabels = labs, asSplits = FALSE)`.
//' @param n_tip Integer; number of tips shared by all trees.
//' @return Numeric vector of length `n*(n-1)/2` containing pairwise MCI
//'   scores in `combn(n, 2)` column-major order (i.e. the data payload of
//'   an R `dist` object).
//' @keywords internal
// [[Rcpp::export]]
NumericVector cpp_mutual_clustering_all_pairs(
    const List& splits_list,
    const int   n_tip,
    const int   n_threads = 1
) {
  const int N = splits_list.size();
  if (N < 2) return NumericVector(0);

  const int n_pairs = N * (N - 1) / 2;

  // Construct SplitList objects on the main thread: R SEXP access is not
  // thread-safe, and the SplitList constructor reads from R RawMatrix objects.
  // SplitList is not move-constructible, so we heap-allocate via unique_ptr.
  std::vector<std::unique_ptr<SplitList>> splits;
  splits.reserve(N);
  for (int k = 0; k < N; ++k) {
    splits.push_back(
      std::make_unique<SplitList>(Rcpp::as<RawMatrix>(splits_list[k]))
    );
  }

  NumericVector result(n_pairs);
  double* const res = result.begin();

  // Iterate over columns of the combn(N,2) lower triangle.
  // Pair (col, row) with col < row maps to dist-vector index:
  //   p = col*(N-1) - col*(col-1)/2 + row - col - 1
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(n_threads)
#endif
  for (int col = 0; col < N - 1; ++col) {
#ifndef _OPENMP
    Rcpp::checkUserInterrupt();
#endif
    for (int row = col + 1; row < N; ++row) {
      const int p = col * (N - 1) - col * (col - 1) / 2 + row - col - 1;
      res[p] = mutual_clustering_score(*splits[col], *splits[row], n_tip);
    }
  }

  return result;
}


// =============================================================================
// InfoRobinsonFoulds (rf_info) — no LAP, thread-safe
// =============================================================================

static double rf_info_score(
    const SplitList& a, const SplitList& b, const int32 n_tips
) {
  const int16 last_bin   = a.n_bins - 1;
  const int16 unset_tips = (n_tips % SL_BIN_SIZE) ?
    SL_BIN_SIZE - n_tips % SL_BIN_SIZE : 0;
  const splitbit unset_mask      = ALL_ONES >> unset_tips;
  const double   lg2_unrooted_n  = lg2_unrooted[n_tips];
  double score = 0;

  splitbit b_complement[SL_MAX_SPLITS][SL_MAX_BINS];
  for (int16 i = 0; i < b.n_splits; ++i) {
    for (int16 bin = 0; bin < last_bin; ++bin) {
      b_complement[i][bin] = ~b.state[i][bin];
    }
    b_complement[i][last_bin] = b.state[i][last_bin] ^ unset_mask;
  }

  for (int16 ai = 0; ai < a.n_splits; ++ai) {
    for (int16 bi = 0; bi < b.n_splits; ++bi) {
      bool all_match = true, all_complement = true;
      for (int16 bin = 0; bin < a.n_bins; ++bin) {
        if (a.state[ai][bin] != b.state[bi][bin]) { all_match = false; break; }
      }
      if (!all_match) {
        for (int16 bin = 0; bin < a.n_bins; ++bin) {
          if (a.state[ai][bin] != b_complement[bi][bin]) {
            all_complement = false; break;
          }
        }
      }
      if (all_match || all_complement) {
        int16 leaves_in_split = 0;
        for (int16 bin = 0; bin < a.n_bins; ++bin) {
          leaves_in_split += count_bits(a.state[ai][bin]);
        }
        score += lg2_unrooted_n
               - lg2_rooted[leaves_in_split]
               - lg2_rooted[n_tips - leaves_in_split];
        break;
      }
    }
  }
  return score;
}

//' @keywords internal
// [[Rcpp::export]]
NumericVector cpp_rf_info_all_pairs(
    const List& splits_list,
    const int   n_tip,
    const int   n_threads = 1
) {
  const int N = splits_list.size();
  if (N < 2) return NumericVector(0);
  const int n_pairs = N * (N - 1) / 2;

  std::vector<std::unique_ptr<SplitList>> splits;
  splits.reserve(N);
  for (int k = 0; k < N; ++k) {
    splits.push_back(
      std::make_unique<SplitList>(Rcpp::as<RawMatrix>(splits_list[k]))
    );
  }

  NumericVector result(n_pairs);
  double* const res = result.begin();

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(n_threads)
#endif
  for (int col = 0; col < N - 1; ++col) {
#ifndef _OPENMP
    Rcpp::checkUserInterrupt();
#endif
    for (int row = col + 1; row < N; ++row) {
      const int p = col * (N - 1) - col * (col - 1) / 2 + row - col - 1;
      res[p] = rf_info_score(*splits[col], *splits[row], n_tip);
    }
  }
  return result;
}


// =============================================================================
// MatchingSplitDistance — LAP-based
// =============================================================================

static double msd_score(
    const SplitList& a, const SplitList& b, const int32 n_tips
) {
  const int16 most_splits = std::max(a.n_splits, b.n_splits);
  if (most_splits == 0) return 0.0;
  const int16 split_diff = most_splits - std::min(a.n_splits, b.n_splits);
  const int16 half_tips  = n_tips / 2;
  const cost  max_score  = BIG / most_splits;

  cost_matrix score(most_splits);

  for (int16 ai = 0; ai < a.n_splits; ++ai) {
    for (int16 bi = 0; bi < b.n_splits; ++bi) {
      splitbit total = 0;
      for (int16 bin = 0; bin < a.n_bins; ++bin) {
        total += count_bits(a.state[ai][bin] ^ b.state[bi][bin]);
      }
      score(ai, bi) = total;
    }
    for (int16 bi = 0; bi < b.n_splits; ++bi) {
      if (score(ai, bi) > half_tips) score(ai, bi) = n_tips - score(ai, bi);
    }
    score.padRowAfterCol(ai, b.n_splits, max_score);
  }
  score.padAfterRow(a.n_splits, max_score);

  std::vector<lap_col> rowsol(most_splits);
  std::vector<lap_row> colsol(most_splits);

  return static_cast<double>(
    lap(most_splits, score, rowsol, colsol, false) - max_score * split_diff
  );
}

//' @keywords internal
// [[Rcpp::export]]
NumericVector cpp_msd_all_pairs(
    const List& splits_list,
    const int   n_tip,
    const int   n_threads = 1
) {
  const int N = splits_list.size();
  if (N < 2) return NumericVector(0);
  const int n_pairs = N * (N - 1) / 2;

  std::vector<std::unique_ptr<SplitList>> splits;
  splits.reserve(N);
  for (int k = 0; k < N; ++k) {
    splits.push_back(
      std::make_unique<SplitList>(Rcpp::as<RawMatrix>(splits_list[k]))
    );
  }

  NumericVector result(n_pairs);
  double* const res = result.begin();

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(n_threads)
#endif
  for (int col = 0; col < N - 1; ++col) {
#ifndef _OPENMP
    Rcpp::checkUserInterrupt();
#endif
    for (int row = col + 1; row < N; ++row) {
      const int p = col * (N - 1) - col * (col - 1) / 2 + row - col - 1;
      res[p] = msd_score(*splits[col], *splits[row], n_tip);
    }
  }
  return result;
}


// =============================================================================
// MatchingSplitInfo — LAP-based
// =============================================================================

static double msi_score(
    const SplitList& a, const SplitList& b, const int32 n_tips
) {
  const int16 most_splits = std::max(a.n_splits, b.n_splits);
  if (most_splits == 0) return 0.0;

  constexpr cost max_score = BIG;
  const double max_possible = lg2_unrooted[n_tips]
    - lg2_rooted[int16((n_tips + 1) / 2)]
    - lg2_rooted[int16(n_tips / 2)];
  const double score_over_possible = static_cast<double>(max_score) / max_possible;
  const double possible_over_score = max_possible / static_cast<double>(max_score);

  cost_matrix score(most_splits);
  splitbit different[SL_MAX_BINS];

  for (int16 ai = 0; ai < a.n_splits; ++ai) {
    for (int16 bi = 0; bi < b.n_splits; ++bi) {
      int16 n_different = 0, n_a_only = 0, n_a_and_b = 0;
      for (int16 bin = 0; bin < a.n_bins; ++bin) {
        different[bin]  = a.state[ai][bin] ^ b.state[bi][bin];
        n_different    += count_bits(different[bin]);
        n_a_only       += count_bits(a.state[ai][bin] &  different[bin]);
        n_a_and_b      += count_bits(a.state[ai][bin] & ~different[bin]);
      }
      const int16 n_same = n_tips - n_different;
      score(ai, bi) = cost(max_score - score_over_possible *
        TreeDist::mmsi_score(n_same, n_a_and_b, n_different, n_a_only));
    }
    score.padRowAfterCol(ai, b.n_splits, max_score);
  }
  score.padAfterRow(a.n_splits, max_score);

  std::vector<lap_col> rowsol(most_splits);
  std::vector<lap_row> colsol(most_splits);

  return static_cast<double>(
    (max_score * most_splits) - lap(most_splits, score, rowsol, colsol, false)
  ) * possible_over_score;
}

//' @keywords internal
// [[Rcpp::export]]
NumericVector cpp_msi_all_pairs(
    const List& splits_list,
    const int   n_tip,
    const int   n_threads = 1
) {
  const int N = splits_list.size();
  if (N < 2) return NumericVector(0);
  const int n_pairs = N * (N - 1) / 2;

  std::vector<std::unique_ptr<SplitList>> splits;
  splits.reserve(N);
  for (int k = 0; k < N; ++k) {
    splits.push_back(
      std::make_unique<SplitList>(Rcpp::as<RawMatrix>(splits_list[k]))
    );
  }

  NumericVector result(n_pairs);
  double* const res = result.begin();

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(n_threads)
#endif
  for (int col = 0; col < N - 1; ++col) {
#ifndef _OPENMP
    Rcpp::checkUserInterrupt();
#endif
    for (int row = col + 1; row < N; ++row) {
      const int p = col * (N - 1) - col * (col - 1) / 2 + row - col - 1;
      res[p] = msi_score(*splits[col], *splits[row], n_tip);
    }
  }
  return result;
}


// =============================================================================
// SharedPhylogeneticInfo — LAP-based
// =============================================================================

static double shared_phylo_score(
    const SplitList& a, const SplitList& b, const int32 n_tips
) {
  const int16 most_splits = std::max(a.n_splits, b.n_splits);
  if (most_splits == 0) return 0.0;

  const int16 overlap_a = int16(n_tips + 1) / 2;
  constexpr cost max_score = BIG;
  const double best_overlap = TreeDist::one_overlap(overlap_a, n_tips / 2, n_tips);
  const double max_possible = lg2_unrooted[n_tips] - best_overlap;
  const double score_over_possible = static_cast<double>(max_score) / max_possible;
  const double possible_over_score = max_possible / static_cast<double>(max_score);

  cost_matrix score(most_splits);

  for (int16 ai = a.n_splits; ai--; ) {
    for (int16 bi = b.n_splits; bi--; ) {
      const double spi = TreeDist::spi_overlap(
        a.state[ai], b.state[bi], n_tips,
        a.in_split[ai], b.in_split[bi], a.n_bins);
      score(ai, bi) = (spi == 0.0) ? max_score
                                   : cost((spi - best_overlap) * score_over_possible);
    }
    score.padRowAfterCol(ai, b.n_splits, max_score);
  }
  score.padAfterRow(a.n_splits, max_score);

  std::vector<lap_col> rowsol(most_splits);
  std::vector<lap_row> colsol(most_splits);

  return static_cast<double>(
    (max_score * most_splits) - lap(most_splits, score, rowsol, colsol, false)
  ) * possible_over_score;
}

//' @keywords internal
// [[Rcpp::export]]
NumericVector cpp_shared_phylo_all_pairs(
    const List& splits_list,
    const int   n_tip,
    const int   n_threads = 1
) {
  const int N = splits_list.size();
  if (N < 2) return NumericVector(0);
  const int n_pairs = N * (N - 1) / 2;

  std::vector<std::unique_ptr<SplitList>> splits;
  splits.reserve(N);
  for (int k = 0; k < N; ++k) {
    splits.push_back(
      std::make_unique<SplitList>(Rcpp::as<RawMatrix>(splits_list[k]))
    );
  }

  NumericVector result(n_pairs);
  double* const res = result.begin();

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(n_threads)
#endif
  for (int col = 0; col < N - 1; ++col) {
#ifndef _OPENMP
    Rcpp::checkUserInterrupt();
#endif
    for (int row = col + 1; row < N; ++row) {
      const int p = col * (N - 1) - col * (col - 1) / 2 + row - col - 1;
      res[p] = shared_phylo_score(*splits[col], *splits[row], n_tip);
    }
  }
  return result;
}


// =============================================================================
// Jaccard / Nye similarity — LAP-based; k and allow_conflict are per-call
// =============================================================================

static double jaccard_score(
    const SplitList& a, const SplitList& b, const int32 n_tips,
    const double exponent, const bool allow_conflict
) {
  const int16 most_splits = std::max(a.n_splits, b.n_splits);
  if (most_splits == 0) return 0.0;

  constexpr cost   max_score  = BIG;
  constexpr double max_scoreL = static_cast<double>(max_score);

  cost_matrix score(most_splits);

  for (int16 ai = 0; ai < a.n_splits; ++ai) {
    const int16 na = a.in_split[ai];
    const int16 nA = n_tips - na;

    for (int16 bi = 0; bi < b.n_splits; ++bi) {
      int16 a_and_b = 0;
      for (int16 bin = 0; bin < a.n_bins; ++bin) {
        a_and_b += count_bits(a.state[ai][bin] & b.state[bi][bin]);
      }
      const int16 nb    = b.in_split[bi];
      const int16 nB    = n_tips - nb;
      const int16 a_and_B = na - a_and_b;
      const int16 A_and_b = nb - a_and_b;
      const int16 A_and_B = nB - a_and_B;

      if (!allow_conflict && !(
            a_and_b == na || a_and_B == na ||
            A_and_b == nA || A_and_B == nA)) {
        score(ai, bi) = max_score;
      } else {
        const int16 A_or_b = n_tips - a_and_B;
        const int16 a_or_B = n_tips - A_and_b;
        const int16 a_or_b = n_tips - A_and_B;
        const int16 A_or_B = n_tips - a_and_b;
        const double ars_ab = static_cast<double>(a_and_b) / static_cast<double>(a_or_b);
        const double ars_Ab = static_cast<double>(A_and_b) / static_cast<double>(A_or_b);
        const double ars_aB = static_cast<double>(a_and_B) / static_cast<double>(a_or_B);
        const double ars_AB = static_cast<double>(A_and_B) / static_cast<double>(A_or_B);
        const double min_both   = (ars_ab < ars_AB) ? ars_ab : ars_AB;
        const double min_either = (ars_aB < ars_Ab) ? ars_aB : ars_Ab;
        const double best = (min_both > min_either) ? min_both : min_either;

        if (exponent == 1.0) {
          score(ai, bi) = cost(max_scoreL - max_scoreL * best);
        } else if (std::isinf(exponent)) {
          score(ai, bi) = cost((best == 1.0) ? 0.0 : max_scoreL);
        } else {
          score(ai, bi) = cost(max_scoreL - max_scoreL * std::pow(best, exponent));
        }
      }
    }
    score.padRowAfterCol(ai, b.n_splits, max_score);
  }
  score.padAfterRow(a.n_splits, max_score);

  std::vector<lap_col> rowsol(most_splits);
  std::vector<lap_row> colsol(most_splits);

  return static_cast<double>(
    (max_score * most_splits) - lap(most_splits, score, rowsol, colsol, false)
  ) / max_scoreL;
}

//' @keywords internal
// [[Rcpp::export]]
NumericVector cpp_jaccard_all_pairs(
    const List&   splits_list,
    const int     n_tip,
    const double  k             = 1.0,
    const bool    allow_conflict = true,
    const int     n_threads     = 1
) {
  const int N = splits_list.size();
  if (N < 2) return NumericVector(0);
  const int n_pairs = N * (N - 1) / 2;

  std::vector<std::unique_ptr<SplitList>> splits;
  splits.reserve(N);
  for (int i = 0; i < N; ++i) {
    splits.push_back(
      std::make_unique<SplitList>(Rcpp::as<RawMatrix>(splits_list[i]))
    );
  }

  NumericVector result(n_pairs);
  double* const res = result.begin();

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(n_threads)
#endif
  for (int col = 0; col < N - 1; ++col) {
#ifndef _OPENMP
    Rcpp::checkUserInterrupt();
#endif
    for (int row = col + 1; row < N; ++row) {
      const int p = col * (N - 1) - col * (col - 1) / 2 + row - col - 1;
      res[p] = jaccard_score(*splits[col], *splits[row], n_tip, k, allow_conflict);
    }
  }
  return result;
}
