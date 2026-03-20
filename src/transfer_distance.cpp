/* transfer_distance.cpp
 *
 * Per-pair and batch computation of the transfer dissimilarity between
 * phylogenetic trees (Takazawa et al. 2026; Lemoine et al. 2018).
 *
 * The transfer distance between two bipartitions is the minimum number of
 * taxa that must be moved to transform one bipartition into the other:
 *   delta(b, b*) = min(Hamming(b, b*), n - Hamming(b, b*))
 *
 * The transfer dissimilarity between two trees sums, for each branch in
 * each tree, the minimum transfer distance to any branch in the other tree.
 */

#include <TreeTools/SplitList.h>
#include <Rcpp/Lightest>
#include <algorithm>
#include <cstdint>
#include <memory>
#include <numeric>
#include <vector>

using Rcpp::IntegerVector;
using Rcpp::List;
using Rcpp::LogicalVector;
using Rcpp::Named;
using Rcpp::NumericMatrix;
using Rcpp::NumericVector;
using Rcpp::RawMatrix;
using TreeTools::count_bits;
using TreeTools::SplitList;

using int16 = int_fast16_t;
using int32 = int_fast32_t;


// ============================================================================
// Per-pair transfer distance computation
// ============================================================================

// Compute transfer distance between two bipartitions (split bitsets).
// Returns min(hamming, n_tips - hamming).
static inline int transfer_dist_splits(
    const splitbit* a, const splitbit* b, int n_bins, int n_tips
) {
  int hamming = 0;
  for (int bin = 0; bin < n_bins; ++bin) {
    hamming += count_bits(a[bin] ^ b[bin]);
  }
  return std::min(hamming, n_tips - hamming);
}


// For each split in `from`, find the minimum transfer distance to any split
// in `to`.  Returns a vector of length from.n_splits with the min distances,
// and a matching vector with the index of the best match in `to` (-1 if no
// improvement over sentinel).
static void find_min_transfer(
    const SplitList& from,
    const SplitList& to,
    int n_tips,
    int n_bins,
    std::vector<int>& min_dist,     // out: [from.n_splits]
    std::vector<int>& best_match    // out: [from.n_splits]
) {
  for (int16 i = 0; i < from.n_splits; ++i) {
    // Sentinel distance: depth(b) - 1
    int pc = 0;
    for (int bin = 0; bin < n_bins; ++bin) {
      pc += count_bits(from.state[i][bin]);
    }
    int depth = std::min(pc, n_tips - pc);
    int sentinel = depth - 1;
    if (sentinel <= 0) {
      min_dist[i] = 0;
      best_match[i] = -1;
      continue;
    }

    int best_d = sentinel;
    int best_j = -1;

    for (int16 j = 0; j < to.n_splits; ++j) {
      int d = transfer_dist_splits(from.state[i], to.state[j], n_bins, n_tips);
      if (d < best_d) {
        best_d = d;
        best_j = j;
        if (d == 0) break; // can't do better
      }
    }

    min_dist[i] = best_d;
    best_match[i] = best_j;
  }
}


// Accumulate the transfer dissimilarity contribution from one direction.
// For scaled: sum of min(min_dist[i] / (depth[i] - 1), 1)
// For unscaled: sum of min(min_dist[i], depth[i] - 1)
static double accumulate_transfer(
    const SplitList& sl,
    const std::vector<int>& min_dist,
    int n_tips, int n_bins,
    bool scale
) {
  double total = 0.0;
  for (int16 i = 0; i < sl.n_splits; ++i) {
    int pc = 0;
    for (int bin = 0; bin < n_bins; ++bin) {
      pc += count_bits(sl.state[i][bin]);
    }
    int depth = std::min(pc, n_tips - pc);
    int p_minus_1 = depth - 1;
    if (p_minus_1 <= 0) continue;

    int d = min_dist[i];
    if (scale) {
      double contrib = static_cast<double>(d) / p_minus_1;
      total += (contrib < 1.0) ? contrib : 1.0;
    } else {
      total += (d < p_minus_1) ? d : p_minus_1;
    }
  }
  return total;
}


// ============================================================================
// Exported: per-pair transfer distance (returns both scaled and unscaled)
// ============================================================================

//' Per-pair transfer dissimilarity
//'
//' @param x,y Raw matrices representing splits (from as.Splits()).
//' @param nTip Integer: number of tips.
//'
//' @return A list with components:
//'   - score_scaled: scaled transfer dissimilarity (double)
//'   - score_unscaled: unscaled transfer dissimilarity (double)
//'   - matching_xy: integer vector, best match in y for each split in x (1-based, NA if sentinel)
//'   - matching_yx: integer vector, best match in x for each split in y (1-based, NA if sentinel)
//' @keywords internal
// [[Rcpp::export]]
List cpp_transfer_dist(
    const RawMatrix& x, const RawMatrix& y,
    const IntegerVector& nTip
) {
  if (x.cols() != y.cols()) {
    Rcpp::stop("Input splits must address same number of tips.");
  }

  const int n_tip = nTip[0];
  SplitList sl_x(x);
  SplitList sl_y(y);
  const int n_bins = sl_x.n_bins;

  // Direction x → y
  std::vector<int> min_dist_xy(sl_x.n_splits);
  std::vector<int> match_xy(sl_x.n_splits);
  find_min_transfer(sl_x, sl_y, n_tip, n_bins, min_dist_xy, match_xy);

  // Direction y → x
  std::vector<int> min_dist_yx(sl_y.n_splits);
  std::vector<int> match_yx(sl_y.n_splits);
  find_min_transfer(sl_y, sl_x, n_tip, n_bins, min_dist_yx, match_yx);

  // Accumulate both directions for scaled and unscaled
  double score_scaled =
      accumulate_transfer(sl_x, min_dist_xy, n_tip, n_bins, true) +
      accumulate_transfer(sl_y, min_dist_yx, n_tip, n_bins, true);
  double score_unscaled =
      accumulate_transfer(sl_x, min_dist_xy, n_tip, n_bins, false) +
      accumulate_transfer(sl_y, min_dist_yx, n_tip, n_bins, false);

  // Build R matching vectors (1-based, NA for sentinel)
  IntegerVector r_match_xy(sl_x.n_splits);
  for (int16 i = 0; i < sl_x.n_splits; ++i) {
    r_match_xy[i] = (match_xy[i] >= 0) ? match_xy[i] + 1 : NA_INTEGER;
  }
  IntegerVector r_match_yx(sl_y.n_splits);
  for (int16 i = 0; i < sl_y.n_splits; ++i) {
    r_match_yx[i] = (match_yx[i] >= 0) ? match_yx[i] + 1 : NA_INTEGER;
  }

  return List::create(
    Named("score_scaled") = score_scaled,
    Named("score_unscaled") = score_unscaled,
    Named("matching_xy") = r_match_xy,
    Named("matching_yx") = r_match_yx
  );
}


// ============================================================================
// Exported: per-pair score-only (for GeneralizedRF compatibility)
// ============================================================================

// Returns a list with "score" and "matching" to be compatible with
// GeneralizedRF / CalculateTreeDistance dispatch.

//' @keywords internal
// [[Rcpp::export]]
List cpp_transfer_dist_scored(
    const RawMatrix& x, const RawMatrix& y,
    const IntegerVector& nTip, bool scale
) {
  if (x.cols() != y.cols()) {
    Rcpp::stop("Input splits must address same number of tips.");
  }

  const int n_tip = nTip[0];
  SplitList sl_x(x);
  SplitList sl_y(y);
  const int n_bins = sl_x.n_bins;

  std::vector<int> min_dist_xy(sl_x.n_splits);
  std::vector<int> match_xy(sl_x.n_splits);
  find_min_transfer(sl_x, sl_y, n_tip, n_bins, min_dist_xy, match_xy);

  std::vector<int> min_dist_yx(sl_y.n_splits);
  std::vector<int> match_yx(sl_y.n_splits);
  find_min_transfer(sl_y, sl_x, n_tip, n_bins, min_dist_yx, match_yx);

  double score = scale
    ? (accumulate_transfer(sl_x, min_dist_xy, n_tip, n_bins, true) +
       accumulate_transfer(sl_y, min_dist_yx, n_tip, n_bins, true))
    : (accumulate_transfer(sl_x, min_dist_xy, n_tip, n_bins, false) +
       accumulate_transfer(sl_y, min_dist_yx, n_tip, n_bins, false));

  // Build matching vector (1-based) — just the x→y direction
  const int n_out = std::max(static_cast<int>(sl_x.n_splits),
                             static_cast<int>(sl_y.n_splits));
  IntegerVector matching(n_out, NA_INTEGER);
  for (int16 i = 0; i < sl_x.n_splits; ++i) {
    matching[i] = (match_xy[i] >= 0) ? match_xy[i] + 1 : NA_INTEGER;
  }

  return List::create(Named("score") = score, Named("matching") = matching);
}


// ============================================================================
// Batch: all-pairs transfer dissimilarity (OpenMP parallel)
// ============================================================================

// Static helper for score-only per-pair computation.
static double transfer_dissimilarity_score(
    const SplitList& sl_x, const SplitList& sl_y,
    int n_tip, int n_bins, bool scale
) {
  // Avoid allocation when trees have 0 splits (star trees)
  if (sl_x.n_splits == 0 && sl_y.n_splits == 0) return 0.0;

  std::vector<int> min_dist_xy(sl_x.n_splits);
  std::vector<int> dummy_match_xy(sl_x.n_splits);
  find_min_transfer(sl_x, sl_y, n_tip, n_bins, min_dist_xy, dummy_match_xy);

  std::vector<int> min_dist_yx(sl_y.n_splits);
  std::vector<int> dummy_match_yx(sl_y.n_splits);
  find_min_transfer(sl_y, sl_x, n_tip, n_bins, min_dist_yx, dummy_match_yx);

  return accumulate_transfer(sl_x, min_dist_xy, n_tip, n_bins, scale) +
         accumulate_transfer(sl_y, min_dist_yx, n_tip, n_bins, scale);
}


//' All-pairs transfer dissimilarity (OpenMP)
//'
//' @param splits_list List of raw matrices (one per tree).
//' @param n_tip Number of tips.
//' @param scale Logical: use scaled transfer dissimilarity?
//' @param n_threads Number of OpenMP threads.
//'
//' @return Numeric vector of length choose(N,2) in dist order.
//' @keywords internal
// [[Rcpp::export]]
NumericVector cpp_transfer_dist_all_pairs(
    const List& splits_list,
    int n_tip, bool scale, int n_threads
) {
  const int N = splits_list.size();
  const int n_pairs = N * (N - 1) / 2;

  // Pre-convert all trees to SplitLists
  std::vector<std::unique_ptr<SplitList>> sls(N);
  int n_bins = 0;
  for (int i = 0; i < N; ++i) {
    sls[i] = std::make_unique<SplitList>(Rcpp::as<RawMatrix>(splits_list[i]));
    if (i == 0) n_bins = sls[0]->n_bins;
  }

  NumericVector result(n_pairs);

  // R dist format: column-major lower triangle
  // pos k corresponds to (row, col) where row > col, iterating col first:
  // (1,0), (2,0), ..., (N-1,0), (2,1), (3,1), ..., (N-1,1), ...
  // k = col * (N - 1) - col * (col - 1) / 2 + (row - col - 1)  (0-based)
  #ifdef _OPENMP
  #pragma omp parallel for schedule(dynamic) num_threads(n_threads)
  #endif
  for (int k = 0; k < n_pairs; ++k) {
    // Invert k to (col, row) in R dist column-major order
    // col = largest c such that c*(2N-c-1)/2 <= k
    int col = 0;
    int remaining = k;
    while (remaining >= N - 1 - col) {
      remaining -= (N - 1 - col);
      ++col;
    }
    int row = col + 1 + remaining;

    result[k] = transfer_dissimilarity_score(
        *sls[col], *sls[row], n_tip, n_bins, scale);
  }

  return result;
}


//' Cross-pairs transfer dissimilarity (OpenMP)
//'
//' @param splits_a,splits_b Lists of raw matrices.
//' @param n_tip Number of tips.
//' @param scale Logical: use scaled transfer dissimilarity?
//' @param n_threads Number of OpenMP threads.
//'
//' @return Numeric matrix of dimension nA x nB.
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix cpp_transfer_dist_cross_pairs(
    const List& splits_a, const List& splits_b,
    int n_tip, bool scale, int n_threads
) {
  const int nA = splits_a.size();
  const int nB = splits_b.size();

  std::vector<std::unique_ptr<SplitList>> sls_a(nA);
  std::vector<std::unique_ptr<SplitList>> sls_b(nB);
  int n_bins = 0;

  for (int i = 0; i < nA; ++i) {
    sls_a[i] = std::make_unique<SplitList>(Rcpp::as<RawMatrix>(splits_a[i]));
    if (i == 0) n_bins = sls_a[0]->n_bins;
  }
  for (int j = 0; j < nB; ++j) {
    sls_b[j] = std::make_unique<SplitList>(Rcpp::as<RawMatrix>(splits_b[j]));
  }

  NumericMatrix result(nA, nB);
  const int total = nA * nB;

  #ifdef _OPENMP
  #pragma omp parallel for schedule(dynamic) num_threads(n_threads)
  #endif
  for (int k = 0; k < total; ++k) {
    int i = k / nB;
    int j = k % nB;
    result(i, j) = transfer_dissimilarity_score(
        *sls_a[i], *sls_b[j], n_tip, n_bins, scale);
  }

  return result;
}
