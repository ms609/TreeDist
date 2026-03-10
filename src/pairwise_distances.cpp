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
