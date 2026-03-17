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
#include <numeric>
#include <vector>
#include <Rcpp/Lightest>

#include "tree_distances.h"

using Rcpp::List;
using Rcpp::NumericVector;
using Rcpp::RawMatrix;
using TreeTools::SplitList;
using TreeTools::count_bits;


// ---------------------------------------------------------------------------
// MatchScratch — reusable buffers for exact-match detection.
// Allocate one per thread and pass into find_exact_matches() to avoid
// per-pair heap allocation.  Vectors grow lazily and are never shrunk.
// ---------------------------------------------------------------------------
struct MatchScratch {
  std::vector<splitbit> a_canon;
  std::vector<splitbit> b_canon;
  std::vector<int16>    a_order;
  std::vector<int16>    b_order;
  std::vector<int16>    a_match;
  std::vector<int16>    b_match;
};

// ---------------------------------------------------------------------------
// O(n log n) exact-match detection for identical bipartitions.
//
// Two splits represent the same bipartition when one equals the other or its
// complement.  We canonicalise each split so that bit 0 (tip 0) is always set;
// then identical bipartitions have identical canonical forms.
//
// Procedure:
//   1. Compute canonical form for every split in a and b  → O(n × n_bins)
//   2. Sort index arrays by canonical form                → O(n log n)
//   3. Merge-scan to find matching pairs                  → O(n)
//
// Returns the number of exact matches found.
// Writes results into scratch.a_match / scratch.b_match:
//   a_match[ai] = bi+1 if split ai matched split bi, else 0.
//   b_match[bi] = ai+1 if split bi matched split ai, else 0.
// ---------------------------------------------------------------------------
static int16 find_exact_matches(
    const SplitList& a, const SplitList& b,
    const int32 n_tips,
    MatchScratch& scratch
) {
  const int16 n_bins   = a.n_bins;
  const int16 last_bin = n_bins - 1;
  const splitbit last_mask = (n_tips % SL_BIN_SIZE == 0)
    ? ~splitbit(0)
    : (splitbit(1) << (n_tips % SL_BIN_SIZE)) - 1;

  const int16 a_n = a.n_splits;
  const int16 b_n = b.n_splits;

  // Ensure buffers are large enough (grow lazily, never shrink)
  const size_t a_canon_sz = static_cast<size_t>(a_n) * n_bins;
  const size_t b_canon_sz = static_cast<size_t>(b_n) * n_bins;
  if (scratch.a_canon.size() < a_canon_sz) scratch.a_canon.resize(a_canon_sz);
  if (scratch.b_canon.size() < b_canon_sz) scratch.b_canon.resize(b_canon_sz);
  if (scratch.a_order.size() < static_cast<size_t>(a_n)) scratch.a_order.resize(a_n);
  if (scratch.b_order.size() < static_cast<size_t>(b_n)) scratch.b_order.resize(b_n);
  if (scratch.a_match.size() < static_cast<size_t>(a_n)) scratch.a_match.resize(a_n);
  if (scratch.b_match.size() < static_cast<size_t>(b_n)) scratch.b_match.resize(b_n);

  int16* a_match = scratch.a_match.data();
  int16* b_match = scratch.b_match.data();
  std::fill(a_match, a_match + a_n, int16(0));
  std::fill(b_match, b_match + b_n, int16(0));

  if (a_n == 0 || b_n == 0) return 0;

  splitbit* a_canon = scratch.a_canon.data();
  splitbit* b_canon = scratch.b_canon.data();

  // --- 1. Compute canonical forms into flat buffers ---
  for (int16 i = 0; i < a_n; ++i) {
    const bool flip = !(a.state[i][0] & 1);
    for (int16 bin = 0; bin < n_bins; ++bin) {
      splitbit val = flip ? ~a.state[i][bin] : a.state[i][bin];
      if (bin == last_bin) val &= last_mask;
      a_canon[i * n_bins + bin] = val;
    }
  }
  for (int16 i = 0; i < b_n; ++i) {
    const bool flip = !(b.state[i][0] & 1);
    for (int16 bin = 0; bin < n_bins; ++bin) {
      splitbit val = flip ? ~b.state[i][bin] : b.state[i][bin];
      if (bin == last_bin) val &= last_mask;
      b_canon[i * n_bins + bin] = val;
    }
  }

  // --- 2. Sort index arrays by canonical form ---
  auto canon_less = [&](const splitbit* canon, int16 i, int16 j) {
    for (int16 bin = 0; bin < n_bins; ++bin) {
      const splitbit vi = canon[i * n_bins + bin];
      const splitbit vj = canon[j * n_bins + bin];
      if (vi < vj) return true;
      if (vi > vj) return false;
    }
    return false; // #nocov
  };

  int16* a_order = scratch.a_order.data();
  int16* b_order = scratch.b_order.data();
  std::iota(a_order, a_order + a_n, int16(0));
  std::iota(b_order, b_order + b_n, int16(0));

  std::sort(a_order, a_order + a_n,
            [&](int16 i, int16 j) {
              return canon_less(a_canon, i, j);
            });
  std::sort(b_order, b_order + b_n,
            [&](int16 i, int16 j) {
              return canon_less(b_canon, i, j);
            });

  // --- 3. Merge-scan to find matches ---
  int16 exact_n = 0;
  int16 ai_pos = 0, bi_pos = 0;
  while (ai_pos < a_n && bi_pos < b_n) {
    const int16 ai = a_order[ai_pos];
    const int16 bi = b_order[bi_pos];

    int cmp = 0;
    for (int16 bin = 0; bin < n_bins; ++bin) {
      const splitbit va = a_canon[ai * n_bins + bin];
      const splitbit vb = b_canon[bi * n_bins + bin];
      if (va < vb) { cmp = -1; break; }
      if (va > vb) { cmp =  1; break; }
    }

    if (cmp < 0) {
      ++ai_pos;
    } else if (cmp > 0) {
      ++bi_pos;
    } else {
      a_match[ai] = bi + 1;
      b_match[bi] = ai + 1;
      ++exact_n;
      ++ai_pos;
      ++bi_pos;
    }
  }

  return exact_n;
}


// Score-only version of mutual_clustering(). Thread-safe: uses only local
// storage and read-only globals (lg2 table, lookup tables).
// Passes allow_interrupt = false to lap() so it is safe to call from an
// OpenMP parallel region.
static double mutual_clustering_score(
    const SplitList& a, const SplitList& b, const int32 n_tips,
    LapScratch& scratch, MatchScratch& mscratch
) {
  if (a.n_splits == 0 || b.n_splits == 0 || n_tips == 0) return 0.0;

  const int16 most_splits = std::max(a.n_splits, b.n_splits);
  const double n_tips_rcp = 1.0 / static_cast<double>(n_tips);

  constexpr cost max_score  = BIG;
  constexpr double over_max = 1.0 / static_cast<double>(BIG);
  const double max_over_tips = static_cast<double>(BIG) * n_tips_rcp;
  const double lg2_n = lg2[n_tips];

  // --- Phase 1: O(n log n) exact-match detection ---
  const int16 exact_n = find_exact_matches(a, b, n_tips, mscratch);
  const int16* a_match = mscratch.a_match.data();
  const int16* b_match = mscratch.b_match.data();

  // Accumulate exact-match score
  double exact_score = 0.0;
  for (int16 ai = 0; ai < a.n_splits; ++ai) {
    if (a_match[ai]) {
      const int16 na = a.in_split[ai];
      const int16 nA = n_tips - na;
      exact_score += TreeDist::ic_matching(na, nA, n_tips);
    }
  }

  // Early exit when everything matched exactly
  if (exact_n == b.n_splits || exact_n == a.n_splits) {
    return exact_score * n_tips_rcp;
  }

  // --- Phase 2: fill cost matrix for unmatched splits only (O(k²)) ---
  const int16 lap_n = most_splits - exact_n;

  // Build index maps for unmatched splits
  std::vector<int16> a_unmatch, b_unmatch;
  a_unmatch.reserve(lap_n);
  b_unmatch.reserve(lap_n);
  for (int16 ai = 0; ai < a.n_splits; ++ai) {
    if (!a_match[ai]) a_unmatch.push_back(ai);
  }
  for (int16 bi = 0; bi < b.n_splits; ++bi) {
    if (!b_match[bi]) b_unmatch.push_back(bi);
  }

  scratch.score_pool.resize(lap_n);
  cost_matrix& score = scratch.score_pool;

  const int16 a_unmatched_n = static_cast<int16>(a_unmatch.size());
  const int16 b_unmatched_n = static_cast<int16>(b_unmatch.size());

  for (int16 a_pos = 0; a_pos < a_unmatched_n; ++a_pos) {
    const int16 ai   = a_unmatch[a_pos];
    const int16 na   = a.in_split[ai];
    const int16 nA   = n_tips - na;
    const auto* a_row = a.state[ai];

    const double offset_a = lg2_n - lg2[na];
    const double offset_A = lg2_n - lg2[nA];

    for (int16 b_pos = 0; b_pos < b_unmatched_n; ++b_pos) {
      const int16 bi   = b_unmatch[b_pos];
      const auto* b_row = b.state[bi];
      int16 a_and_b = 0;
      for (int16 bin = 0; bin < a.n_bins; ++bin) {
        a_and_b += count_bits(a_row[bin] & b_row[bin]);
      }

      const int16 nb    = b.in_split[bi];
      const int16 nB    = n_tips - nb;
      const int16 a_and_B = na - a_and_b;
      const int16 A_and_b = nb - a_and_b;
      const int16 A_and_B = nA - A_and_b;

      if (a_and_b == A_and_b && a_and_b == a_and_B && a_and_b == A_and_B) {
        score(a_pos, b_pos) = max_score;
      } else {
        const double lg2_nb = lg2[nb];
        const double lg2_nB = lg2[nB];
        const double ic_sum =
          a_and_b * (lg2[a_and_b] + offset_a - lg2_nb) +
          a_and_B * (lg2[a_and_B] + offset_a - lg2_nB) +
          A_and_b * (lg2[A_and_b] + offset_A - lg2_nb) +
          A_and_B * (lg2[A_and_B] + offset_A - lg2_nB);
        score(a_pos, b_pos) = max_score - static_cast<cost>(ic_sum * max_over_tips);
      }
    }
    // Pad extra columns for asymmetric split counts
    if (b_unmatched_n < lap_n) {
      score.padRowAfterCol(a_pos, b_unmatched_n, max_score);
    }
  }
  // Pad extra rows for asymmetric split counts
  if (a_unmatched_n < lap_n) {
    score.padAfterRow(a_unmatched_n, max_score);
  }

  // --- Phase 3: solve LAP on the reduced matrix ---
  scratch.ensure(lap_n);
  auto& rowsol = scratch.rowsol;
  auto& colsol = scratch.colsol;

  const double lap_score =
    static_cast<double>((max_score * lap_n) -
                        lap(lap_n, score, rowsol, colsol, false, scratch)) * over_max;
  return lap_score + exact_score * n_tips_rcp;
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

  // One scratch set per thread — grown lazily on first use, never freed between
  // pairs.  Indexed by omp_get_thread_num() (always 0 in the serial path).
  const int n_scratch = std::max(1, n_threads);
  std::vector<LapScratch> scratches(n_scratch);
  std::vector<MatchScratch> mscratches(n_scratch);

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
#ifdef _OPENMP
    const int tid = omp_get_thread_num();
#else
    const int tid = 0;
#endif
    LapScratch& scratch = scratches[tid];
    MatchScratch& mscratch = mscratches[tid];
    for (int row = col + 1; row < N; ++row) {
      const int p = col * (N - 1) - col * (col - 1) / 2 + row - col - 1;
      res[p] = mutual_clustering_score(*splits[col], *splits[row], n_tip, scratch, mscratch);
    }
  }

  return result;
}


// =============================================================================
// InfoRobinsonFoulds (rf_info) — no LAP, thread-safe
// =============================================================================

static double rf_info_score(
    const SplitList& a, const SplitList& b, const int32 n_tips,
    MatchScratch& mscratch
) {
  const int16 a_n = a.n_splits;
  const int16 b_n = b.n_splits;
  if (a_n == 0 || b_n == 0) return 0;

  // Use sort+merge to find exact matches in O(n log n)
  const int16 exact_n = find_exact_matches(a, b, n_tips, mscratch);
  if (exact_n == 0) return 0;

  // Sum info contribution for each matched split in a
  const int16* a_match = mscratch.a_match.data();
  const double lg2_unrooted_n = lg2_unrooted[n_tips];
  double score = 0;
  for (int16 ai = 0; ai < a_n; ++ai) {
    if (a_match[ai] == 0) continue;
    int16 leaves_in_split = 0;
    for (int16 bin = 0; bin < a.n_bins; ++bin) {
      leaves_in_split += count_bits(a.state[ai][bin]);
    }
    score += lg2_unrooted_n
           - lg2_rooted[leaves_in_split]
           - lg2_rooted[n_tips - leaves_in_split];
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

  const int n_scratch = std::max(1, n_threads);
  std::vector<MatchScratch> mscratches(n_scratch);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(n_threads)
#endif
  for (int col = 0; col < N - 1; ++col) {
#ifndef _OPENMP
    Rcpp::checkUserInterrupt();
#endif
#ifdef _OPENMP
    MatchScratch& mscratch = mscratches[omp_get_thread_num()];
#else
    MatchScratch& mscratch = mscratches[0];
#endif
    for (int row = col + 1; row < N; ++row) {
      const int p = col * (N - 1) - col * (col - 1) / 2 + row - col - 1;
      res[p] = rf_info_score(*splits[col], *splits[row], n_tip, mscratch);
    }
  }
  return result;
}


// =============================================================================
// MatchingSplitDistance — LAP-based
// =============================================================================

static double msd_score(
    const SplitList& a, const SplitList& b, const int32 n_tips,
    LapScratch& scratch, MatchScratch& mscratch
) {
  const int16 most_splits = std::max(a.n_splits, b.n_splits);
  if (most_splits == 0) return 0.0;
  const bool  a_has_more  = (a.n_splits > b.n_splits);
  const int16 a_extra     = a_has_more ? most_splits - b.n_splits : 0;
  const int16 b_extra     = a_has_more ? 0 : most_splits - a.n_splits;
  const int16 half_tips   = n_tips / 2;
  const cost  max_score   = BIG / most_splits;

  // --- Phase 1: O(n log n) exact-match detection ---
  const int16 exact_n = find_exact_matches(a, b, n_tips, mscratch);
  const int16* a_match = mscratch.a_match.data();
  const int16* b_match = mscratch.b_match.data();

  if (exact_n == b.n_splits || exact_n == a.n_splits) {
    return 0.0;
  }

  // --- Phase 2: fill cost matrix for unmatched splits only ---
  const int16 lap_n = most_splits - exact_n;

  std::vector<int16> a_unmatch, b_unmatch;
  a_unmatch.reserve(lap_n);
  b_unmatch.reserve(lap_n);
  for (int16 ai = 0; ai < a.n_splits; ++ai) {
    if (!a_match[ai]) a_unmatch.push_back(ai);
  }
  for (int16 bi = 0; bi < b.n_splits; ++bi) {
    if (!b_match[bi]) b_unmatch.push_back(bi);
  }

  scratch.score_pool.resize(lap_n);
  cost_matrix& score = scratch.score_pool;

  const int16 a_unmatched_n = static_cast<int16>(a_unmatch.size());
  const int16 b_unmatched_n = static_cast<int16>(b_unmatch.size());

  for (int16 a_pos = 0; a_pos < a_unmatched_n; ++a_pos) {
    const int16 ai = a_unmatch[a_pos];
    for (int16 b_pos = 0; b_pos < b_unmatched_n; ++b_pos) {
      const int16 bi = b_unmatch[b_pos];
      splitbit total = 0;
      for (int16 bin = 0; bin < a.n_bins; ++bin) {
        total += count_bits(a.state[ai][bin] ^ b.state[bi][bin]);
      }
      score(a_pos, b_pos) = static_cast<cost>(
        total > static_cast<splitbit>(half_tips) ? n_tips - total : total);
    }
    if (b_unmatched_n < lap_n) {
      score.padRowAfterCol(a_pos, b_unmatched_n, max_score);
    }
  }
  if (a_unmatched_n < lap_n) {
    score.padAfterRow(a_unmatched_n, max_score);
  }

  // --- Phase 3: solve LAP ---
  scratch.ensure(lap_n);
  auto& rowsol = scratch.rowsol;
  auto& colsol = scratch.colsol;

  const cost split_diff_cost = max_score * (a_extra + b_extra);
  return static_cast<double>(
    lap(lap_n, score, rowsol, colsol, false, scratch) - split_diff_cost);
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

  const int n_scratch = std::max(1, n_threads);
  std::vector<LapScratch> scratches(n_scratch);
  std::vector<MatchScratch> mscratches(n_scratch);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(n_threads)
#endif
  for (int col = 0; col < N - 1; ++col) {
#ifndef _OPENMP
    Rcpp::checkUserInterrupt();
#endif
#ifdef _OPENMP
    const int tid = omp_get_thread_num();
#else
    const int tid = 0;
#endif
    LapScratch& scratch = scratches[tid];
    MatchScratch& mscratch = mscratches[tid];
    for (int row = col + 1; row < N; ++row) {
      const int p = col * (N - 1) - col * (col - 1) / 2 + row - col - 1;
      res[p] = msd_score(*splits[col], *splits[row], n_tip, scratch, mscratch);
    }
  }
  return result;
}


// =============================================================================
// MatchingSplitInfo — LAP-based
// =============================================================================

// NOTE: msi_score does NOT use exact-match detection.  For the MSI metric,
// matching a split to its identical copy is NOT necessarily globally optimal —
// the LAP can find a better assignment where the identical split pairs with
// a different (compatible, containing) split that has higher mmsi_score.
// This differs from MCI and MSD where exact matches are provably optimal.
static double msi_score(
    const SplitList& a, const SplitList& b, const int32 n_tips,
    LapScratch& scratch
) {
  const int16 most_splits = std::max(a.n_splits, b.n_splits);
  if (most_splits == 0) return 0.0;

  constexpr cost max_score = BIG;
  const double max_possible = lg2_unrooted[n_tips]
    - lg2_rooted[int16((n_tips + 1) / 2)]
    - lg2_rooted[int16(n_tips / 2)];
  const double score_over_possible = static_cast<double>(max_score) / max_possible;
  const double possible_over_score = max_possible / static_cast<double>(max_score);

  scratch.score_pool.resize(most_splits);
  cost_matrix& score = scratch.score_pool;

  for (int16 ai = 0; ai < a.n_splits; ++ai) {
    for (int16 bi = 0; bi < b.n_splits; ++bi) {
      int16 n_a_only = 0, n_a_and_b = 0, n_different = 0;
      splitbit different;
      for (int16 bin = 0; bin < a.n_bins; ++bin) {
        different   = a.state[ai][bin] ^ b.state[bi][bin];
        n_different += count_bits(different);
        n_a_only   += count_bits(a.state[ai][bin] &  different);
        n_a_and_b  += count_bits(a.state[ai][bin] & ~different);
      }
      const int16 n_same = n_tips - n_different;
      score(ai, bi) = cost(max_score - score_over_possible *
        TreeDist::mmsi_score(n_same, n_a_and_b, n_different, n_a_only));
    }
    score.padRowAfterCol(ai, b.n_splits, max_score);
  }
  score.padAfterRow(a.n_splits, max_score);

  scratch.ensure(most_splits);
  auto& rowsol = scratch.rowsol;
  auto& colsol = scratch.colsol;

  return static_cast<double>(
    (max_score * most_splits) - lap(most_splits, score, rowsol, colsol, false, scratch)
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

  const int n_scratch = std::max(1, n_threads);
  std::vector<LapScratch> scratches(n_scratch);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(n_threads)
#endif
  for (int col = 0; col < N - 1; ++col) {
#ifndef _OPENMP
    Rcpp::checkUserInterrupt();
#endif
#ifdef _OPENMP
    LapScratch& scratch = scratches[omp_get_thread_num()];
#else
    LapScratch& scratch = scratches[0];
#endif
    for (int row = col + 1; row < N; ++row) {
      const int p = col * (N - 1) - col * (col - 1) / 2 + row - col - 1;
      res[p] = msi_score(*splits[col], *splits[row], n_tip, scratch);
    }
  }
  return result;
}


// =============================================================================
// SharedPhylogeneticInfo — LAP-based
// =============================================================================

// NOTE: shared_phylo_score does NOT use exact-match detection.  For the SPI
// metric, spi_overlap(A, B) where B contains A can EXCEED spi_overlap(A, A)
// for identical splits, so the LAP can find a better global assignment by NOT
// matching identical splits.  The full LAP is always solved.
static double shared_phylo_score(
    const SplitList& a, const SplitList& b, const int32 n_tips,
    LapScratch& scratch
) {
  const int16 most_splits = std::max(a.n_splits, b.n_splits);
  if (most_splits == 0) return 0.0;

  const int16 overlap_a = int16(n_tips + 1) / 2;
  constexpr cost max_score = BIG;
  const double best_overlap = TreeDist::one_overlap(overlap_a, n_tips / 2, n_tips);
  const double max_possible = lg2_unrooted[n_tips] - best_overlap;
  const double score_over_possible = static_cast<double>(max_score) / max_possible;
  const double possible_over_score = max_possible / static_cast<double>(max_score);

  scratch.score_pool.resize(most_splits);
  cost_matrix& score = scratch.score_pool;

  for (int16 ai = 0; ai < a.n_splits; ++ai) {
    for (int16 bi = 0; bi < b.n_splits; ++bi) {
      const double spi = TreeDist::spi_overlap(
        a.state[ai], b.state[bi], n_tips,
        a.in_split[ai], b.in_split[bi], a.n_bins);
      score(ai, bi) = (spi == 0.0) ? max_score
                                   : cost((spi - best_overlap) * score_over_possible);
    }
    score.padRowAfterCol(ai, b.n_splits, max_score);
  }
  score.padAfterRow(a.n_splits, max_score);

  scratch.ensure(most_splits);
  auto& rowsol = scratch.rowsol;
  auto& colsol = scratch.colsol;

  return static_cast<double>(
    (max_score * most_splits) - lap(most_splits, score, rowsol, colsol, false, scratch)
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

  const int n_scratch = std::max(1, n_threads);
  std::vector<LapScratch> scratches(n_scratch);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(n_threads)
#endif
  for (int col = 0; col < N - 1; ++col) {
#ifndef _OPENMP
    Rcpp::checkUserInterrupt();
#endif
#ifdef _OPENMP
    LapScratch& scratch = scratches[omp_get_thread_num()];
#else
    LapScratch& scratch = scratches[0];
#endif
    for (int row = col + 1; row < N; ++row) {
      const int p = col * (N - 1) - col * (col - 1) / 2 + row - col - 1;
      res[p] = shared_phylo_score(*splits[col], *splits[row], n_tip, scratch);
    }
  }
  return result;
}


// =============================================================================
// Jaccard / Nye similarity — LAP-based; k and allow_conflict are per-call
// =============================================================================

static double jaccard_score(
    const SplitList& a, const SplitList& b, const int32 n_tips,
    const double exponent, const bool allow_conflict,
    LapScratch& scratch, MatchScratch& mscratch
) {
  const int16 most_splits = std::max(a.n_splits, b.n_splits);
  if (most_splits == 0) return 0.0;

  constexpr cost   max_score  = BIG;
  constexpr double max_scoreL = static_cast<double>(max_score);

  // --- Phase 1: O(n log n) exact-match detection ---
  // Only used when allow_conflict=true; otherwise the full LAP may reassign
  // non-matching splits to compatible (non-exact) partners.
  int16 exact_n = 0;
  if (allow_conflict) {
    exact_n = find_exact_matches(a, b, n_tips, mscratch);
  } else {
    // Ensure match arrays are sized and zeroed
    if (mscratch.a_match.size() < static_cast<size_t>(a.n_splits))
      mscratch.a_match.resize(a.n_splits);
    if (mscratch.b_match.size() < static_cast<size_t>(b.n_splits))
      mscratch.b_match.resize(b.n_splits);
    std::fill(mscratch.a_match.data(), mscratch.a_match.data() + a.n_splits, int16(0));
    std::fill(mscratch.b_match.data(), mscratch.b_match.data() + b.n_splits, int16(0));
  }
  const int16* a_match = mscratch.a_match.data();
  const int16* b_match = mscratch.b_match.data();

  if (exact_n == b.n_splits || exact_n == a.n_splits) {
    return static_cast<double>(exact_n);
  }

  // --- Phase 2: fill cost matrix for unmatched splits only ---
  const int16 lap_n = most_splits - exact_n;

  std::vector<int16> a_unmatch, b_unmatch;
  a_unmatch.reserve(lap_n);
  b_unmatch.reserve(lap_n);
  for (int16 ai = 0; ai < a.n_splits; ++ai) {
    if (!a_match[ai]) a_unmatch.push_back(ai);
  }
  for (int16 bi = 0; bi < b.n_splits; ++bi) {
    if (!b_match[bi]) b_unmatch.push_back(bi);
  }

  scratch.score_pool.resize(lap_n);
  cost_matrix& score = scratch.score_pool;

  const int16 a_unmatched_n = static_cast<int16>(a_unmatch.size());
  const int16 b_unmatched_n = static_cast<int16>(b_unmatch.size());

  for (int16 a_pos = 0; a_pos < a_unmatched_n; ++a_pos) {
    const int16 ai = a_unmatch[a_pos];
    const int16 na = a.in_split[ai];
    const int16 nA = n_tips - na;

    for (int16 b_pos = 0; b_pos < b_unmatched_n; ++b_pos) {
      const int16 bi = b_unmatch[b_pos];
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
        score(a_pos, b_pos) = max_score;
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
          score(a_pos, b_pos) = cost(max_scoreL - max_scoreL * best);
        } else if (std::isinf(exponent)) {
          score(a_pos, b_pos) = cost((best == 1.0) ? 0.0 : max_scoreL);
        } else {
          score(a_pos, b_pos) = cost(max_scoreL - max_scoreL * std::pow(best, exponent));
        }
      }
    }
    if (b_unmatched_n < lap_n) {
      score.padRowAfterCol(a_pos, b_unmatched_n, max_score);
    }
  }
  if (a_unmatched_n < lap_n) {
    score.padAfterRow(a_unmatched_n, max_score);
  }

  // --- Phase 3: solve LAP ---
  scratch.ensure(lap_n);
  auto& rowsol = scratch.rowsol;
  auto& colsol = scratch.colsol;

  return static_cast<double>(exact_n) + static_cast<double>(
    (max_score * lap_n) - lap(lap_n, score, rowsol, colsol, false, scratch)
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

  const int n_scratch = std::max(1, n_threads);
  std::vector<LapScratch> scratches(n_scratch);
  std::vector<MatchScratch> mscratches(n_scratch);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(n_threads)
#endif
  for (int col = 0; col < N - 1; ++col) {
#ifndef _OPENMP
    Rcpp::checkUserInterrupt();
#endif
#ifdef _OPENMP
    const int tid = omp_get_thread_num();
#else
    const int tid = 0;
#endif
    LapScratch& scratch = scratches[tid];
    MatchScratch& mscratch = mscratches[tid];
    for (int row = col + 1; row < N; ++row) {
      const int p = col * (N - 1) - col * (col - 1) / 2 + row - col - 1;
      res[p] = jaccard_score(*splits[col], *splits[row], n_tip, k, allow_conflict, scratch, mscratch);
    }
  }
  return result;
}


// =============================================================================
// Cross-pairs (ManyMany) batch functions
// =============================================================================
//
// Each function computes the full nA × nB matrix of pairwise distances between
// two sets of trees.  Like the all_pairs variants, these use OpenMP where
// available and reuse LapScratch / CostMatrix pools across pairs.

static void parse_split_list(const List& list,
                             std::vector<std::unique_ptr<SplitList>>& out) {
  const int n = list.size();
  out.reserve(n);
  for (int k = 0; k < n; ++k) {
    out.push_back(
      std::make_unique<SplitList>(Rcpp::as<RawMatrix>(list[k]))
    );
  }
}

using Rcpp::NumericMatrix;

//' @keywords internal
// [[Rcpp::export]]
NumericMatrix cpp_mutual_clustering_cross_pairs(
    const List& splits_a, const List& splits_b,
    const int n_tip, const int n_threads = 1
) {
  const int nA = splits_a.size();
  const int nB = splits_b.size();
  if (nA == 0 || nB == 0) return NumericMatrix(nA, nB);

  std::vector<std::unique_ptr<SplitList>> sa, sb;
  parse_split_list(splits_a, sa);
  parse_split_list(splits_b, sb);

  NumericMatrix result(nA, nB);
  double* res = result.begin();

  const int n_scratch = std::max(1, n_threads);
  std::vector<LapScratch> scratches(n_scratch);
  std::vector<MatchScratch> mscratches(n_scratch);
  const int total = nA * nB;

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(n_threads)
#endif
  for (int idx = 0; idx < total; ++idx) {
#ifndef _OPENMP
    if ((idx & 0xFF) == 0) Rcpp::checkUserInterrupt();
#endif
#ifdef _OPENMP
    const int tid = omp_get_thread_num();
#else
    const int tid = 0;
#endif
    LapScratch& scratch = scratches[tid];
    MatchScratch& mscratch = mscratches[tid];
    const int i = idx % nA;
    const int j = idx / nA;
    res[idx] = mutual_clustering_score(*sa[i], *sb[j], n_tip, scratch, mscratch);
  }
  return result;
}

//' @keywords internal
// [[Rcpp::export]]
NumericMatrix cpp_rf_info_cross_pairs(
    const List& splits_a, const List& splits_b,
    const int n_tip, const int n_threads = 1
) {
  const int nA = splits_a.size();
  const int nB = splits_b.size();
  if (nA == 0 || nB == 0) return NumericMatrix(nA, nB);

  std::vector<std::unique_ptr<SplitList>> sa, sb;
  parse_split_list(splits_a, sa);
  parse_split_list(splits_b, sb);

  NumericMatrix result(nA, nB);
  double* res = result.begin();
  const int total = nA * nB;

  const int n_scratch = std::max(1, n_threads);
  std::vector<MatchScratch> mscratches(n_scratch);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(n_threads)
#endif
  for (int idx = 0; idx < total; ++idx) {
#ifndef _OPENMP
    if ((idx & 0xFF) == 0) Rcpp::checkUserInterrupt();
#endif
#ifdef _OPENMP
    MatchScratch& mscratch = mscratches[omp_get_thread_num()];
#else
    MatchScratch& mscratch = mscratches[0];
#endif
    const int i = idx % nA;
    const int j = idx / nA;
    res[idx] = rf_info_score(*sa[i], *sb[j], n_tip, mscratch);
  }
  return result;
}

//' @keywords internal
// [[Rcpp::export]]
NumericMatrix cpp_msd_cross_pairs(
    const List& splits_a, const List& splits_b,
    const int n_tip, const int n_threads = 1
) {
  const int nA = splits_a.size();
  const int nB = splits_b.size();
  if (nA == 0 || nB == 0) return NumericMatrix(nA, nB);

  std::vector<std::unique_ptr<SplitList>> sa, sb;
  parse_split_list(splits_a, sa);
  parse_split_list(splits_b, sb);

  NumericMatrix result(nA, nB);
  double* res = result.begin();

  const int n_scratch = std::max(1, n_threads);
  std::vector<LapScratch> scratches(n_scratch);
  std::vector<MatchScratch> mscratches(n_scratch);
  const int total = nA * nB;

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(n_threads)
#endif
  for (int idx = 0; idx < total; ++idx) {
#ifndef _OPENMP
    if ((idx & 0xFF) == 0) Rcpp::checkUserInterrupt();
#endif
#ifdef _OPENMP
    const int tid = omp_get_thread_num();
#else
    const int tid = 0;
#endif
    LapScratch& scratch = scratches[tid];
    MatchScratch& mscratch = mscratches[tid];
    const int i = idx % nA;
    const int j = idx / nA;
    res[idx] = msd_score(*sa[i], *sb[j], n_tip, scratch, mscratch);
  }
  return result;
}

//' @keywords internal
// [[Rcpp::export]]
NumericMatrix cpp_msi_cross_pairs(
    const List& splits_a, const List& splits_b,
    const int n_tip, const int n_threads = 1
) {
  const int nA = splits_a.size();
  const int nB = splits_b.size();
  if (nA == 0 || nB == 0) return NumericMatrix(nA, nB);

  std::vector<std::unique_ptr<SplitList>> sa, sb;
  parse_split_list(splits_a, sa);
  parse_split_list(splits_b, sb);

  NumericMatrix result(nA, nB);
  double* res = result.begin();

  const int n_scratch = std::max(1, n_threads);
  std::vector<LapScratch> scratches(n_scratch);
  const int total = nA * nB;

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(n_threads)
#endif
  for (int idx = 0; idx < total; ++idx) {
#ifndef _OPENMP
    if ((idx & 0xFF) == 0) Rcpp::checkUserInterrupt();
#endif
#ifdef _OPENMP
    LapScratch& scratch = scratches[omp_get_thread_num()];
#else
    LapScratch& scratch = scratches[0];
#endif
    const int i = idx % nA;
    const int j = idx / nA;
    res[idx] = msi_score(*sa[i], *sb[j], n_tip, scratch);
  }
  return result;
}

//' @keywords internal
// [[Rcpp::export]]
NumericMatrix cpp_shared_phylo_cross_pairs(
    const List& splits_a, const List& splits_b,
    const int n_tip, const int n_threads = 1
) {
  const int nA = splits_a.size();
  const int nB = splits_b.size();
  if (nA == 0 || nB == 0) return NumericMatrix(nA, nB);

  std::vector<std::unique_ptr<SplitList>> sa, sb;
  parse_split_list(splits_a, sa);
  parse_split_list(splits_b, sb);

  NumericMatrix result(nA, nB);
  double* res = result.begin();

  const int n_scratch = std::max(1, n_threads);
  std::vector<LapScratch> scratches(n_scratch);
  const int total = nA * nB;

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(n_threads)
#endif
  for (int idx = 0; idx < total; ++idx) {
#ifndef _OPENMP
    if ((idx & 0xFF) == 0) Rcpp::checkUserInterrupt();
#endif
#ifdef _OPENMP
    LapScratch& scratch = scratches[omp_get_thread_num()];
#else
    LapScratch& scratch = scratches[0];
#endif
    const int i = idx % nA;
    const int j = idx / nA;
    res[idx] = shared_phylo_score(*sa[i], *sb[j], n_tip, scratch);
  }
  return result;
}

//' @keywords internal
// [[Rcpp::export]]
NumericMatrix cpp_jaccard_cross_pairs(
    const List& splits_a, const List& splits_b,
    const int n_tip,
    const double k = 1.0,
    const bool allow_conflict = true,
    const int n_threads = 1
) {
  const int nA = splits_a.size();
  const int nB = splits_b.size();
  if (nA == 0 || nB == 0) return NumericMatrix(nA, nB);

  std::vector<std::unique_ptr<SplitList>> sa, sb;
  parse_split_list(splits_a, sa);
  parse_split_list(splits_b, sb);

  NumericMatrix result(nA, nB);
  double* res = result.begin();

  const int n_scratch = std::max(1, n_threads);
  std::vector<LapScratch> scratches(n_scratch);
  std::vector<MatchScratch> mscratches(n_scratch);
  const int total = nA * nB;

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(n_threads)
#endif
  for (int idx = 0; idx < total; ++idx) {
#ifndef _OPENMP
    if ((idx & 0xFF) == 0) Rcpp::checkUserInterrupt();
#endif
#ifdef _OPENMP
    const int tid = omp_get_thread_num();
#else
    const int tid = 0;
#endif
    LapScratch& scratch = scratches[tid];
    MatchScratch& mscratch = mscratches[tid];
    const int i = idx % nA;
    const int j = idx / nA;
    res[idx] = jaccard_score(*sa[i], *sb[j], n_tip, k, allow_conflict, scratch, mscratch);
  }
  return result;
}


// ---------------------------------------------------------------------------
// Per-tree entropy / information batch functions
//
// These replace per-tree R vapply() calls in .FastDistPath() with a single
// C++ call over all trees, eliminating N R->C++ round-trips.
// ---------------------------------------------------------------------------

// ClusteringEntropy: sum of binary_entropy(k/n) over splits, in nats then
// converted to bits (divide by log(2)).
// binary_entropy(p) = -(p log p + (1-p) log(1-p))
// [[Rcpp::export]]
NumericVector cpp_clustering_entropy_batch(
    const List& splits_list,
    const int   n_tip
) {
  const int N = splits_list.size();
  NumericVector result(N);
  if (N == 0 || n_tip <= 0) return result;

  const double invN = 1.0 / static_cast<double>(n_tip);
  constexpr double invLog2 = 1.442695040888963387005;

  for (int i = 0; i < N; ++i) {
    SplitList sl(Rcpp::as<RawMatrix>(splits_list[i]));
    double total = 0.0;
    for (int16 s = 0; s < sl.n_splits; ++s) {
      const int k = sl.in_split[s];
      if (k <= 0 || k >= n_tip) continue;
      const double p = k * invN;
      const double q = 1.0 - p;
      total += -(p * std::log(p) + q * std::log1p(-p));
    }
    result[i] = total * invLog2;
  }
  return result;
}

// SplitwiseInfo: sum of (Log2Unrooted(n) - Log2Rooted(k) - Log2Rooted(n-k))
// over splits.
// Log2Rooted(m) = log2((2m-3)!!) = sum_{j=1,3,...,2m-3} log2(j)
// Log2Unrooted(n) = Log2Rooted(n) - log2(2n-3)  [= log2((2n-5)!!)]
//
// We precompute a table of Log2Rooted values up to n_tip.
// [[Rcpp::export]]
NumericVector cpp_splitwise_info_batch(
    const List& splits_list,
    const int   n_tip
) {
  const int N = splits_list.size();
  NumericVector result(N);
  if (N == 0 || n_tip < 4) return result;

  // Precompute Log2Rooted table: log2((2m-3)!!) for m = 0..n_tip
  // Log2Rooted(0) = Log2Rooted(1) = 0 (convention)
  // Log2Rooted(2) = log2(1) = 0
  // Log2Rooted(m) = Log2Rooted(m-1) + log2(2m-3)
  std::vector<double> l2r(n_tip + 1, 0.0);
  for (int m = 3; m <= n_tip; ++m) {
    l2r[m] = l2r[m - 1] + std::log2(2 * m - 3);
  }
  // Log2Unrooted(n) = log2((2n-5)!!) = Log2Rooted(n) - log2(2n-3)
  const double l2u_n = l2r[n_tip] - std::log2(2 * n_tip - 3);

  for (int i = 0; i < N; ++i) {
    SplitList sl(Rcpp::as<RawMatrix>(splits_list[i]));
    double total = 0.0;
    for (int16 s = 0; s < sl.n_splits; ++s) {
      const int k = sl.in_split[s];
      if (k < 2 || (n_tip - k) < 2) continue;
      total += l2u_n - l2r[k] - l2r[n_tip - k];
    }
    result[i] = total;
  }
  return result;
}
