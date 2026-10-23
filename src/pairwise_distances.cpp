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
#include <TreeDist/popcount.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <cstring>
#include <numeric>
#include <vector>
#include <Rcpp/Lightest>

#include "tree_distances.h"

using Rcpp::List;
using Rcpp::NumericVector;
using Rcpp::RawMatrix;
using TreeTools::SplitList;


// ---------------------------------------------------------------------------
// MatchScratch — reusable buffers for exact-match detection.
// Allocate one per thread and pass into find_exact_matches() to avoid
// per-pair heap allocation.  Vectors grow lazily and are never shrunk.
// ---------------------------------------------------------------------------
struct MatchScratch {
  std::vector<split_int>  a_match;
  std::vector<split_int>  b_match;
};

// Reusable buffers for shared_phylo_pinned(); one per thread.
struct PinScratch {
  std::vector<split_int> pin_a, pin_b, free_a, free_b;
  std::vector<cost> pin_d, pin_lo, pin_u, u;
  std::vector<char> b_pinned;
  std::vector<split_int> cls;          // size group of each pinned pair
  std::vector<split_int> group_of;     // group index by split size, or -1
  std::vector<split_int> group_size, group_start, by_u, by_w;
};

// ---------------------------------------------------------------------------
// CanonSplits — per-tree canonical split forms, held in sorted order.
//
// Two splits represent the same bipartition when one equals the other or its
// complement.  Canonicalising each split so that bit 0 (tip 0) is always set
// makes identical bipartitions compare equal.
//
// The canonical form and its sort order depend on one tree only, so they are
// computed once per tree (build_canon_all) rather than once per pair: a set of
// N trees needs N such builds, not the N(N-1) implied by building them inside
// find_exact_matches().
//
// `canon` stores the canonical forms already permuted into sorted order, so the
// merge scan reads two contiguous arrays instead of gathering through an index
// vector.  `idx[pos]` recovers the original split number at sorted position
// `pos`.
//
// Memory: n_splits × n_bins × sizeof(splitbit) per tree — the same footprint as
// the SplitList `state` the caller already holds, so retaining these for a whole
// tree set doubles split storage.
// ---------------------------------------------------------------------------
struct CanonSplits {
  std::vector<splitbit>  canon;  // n_splits × n_bins, in sorted order
  std::vector<split_int> idx;    // idx[pos] = original split index
};

static void build_canon(const SplitList& s, const int32 n_tips,
                        CanonSplits& out) {
  const split_int n        = s.n_splits;
  const split_int n_bins   = s.n_bins;
  const split_int last_bin = n_bins - 1;
  const splitbit last_mask = (n_tips % SL_BIN_SIZE == 0)
    ? ~splitbit(0)
    : (splitbit(1) << (n_tips % SL_BIN_SIZE)) - 1;

  out.canon.resize(static_cast<size_t>(n) * n_bins);
  out.idx.resize(n);
  if (n == 0) return;

  // --- 1. Canonical form, indexed by original split number ---
  std::vector<splitbit> tmp(static_cast<size_t>(n) * n_bins);
  for (split_int i = 0; i < n; ++i) {
    const bool flip = !(s.state[i][0] & 1);
    for (split_int bin = 0; bin < n_bins; ++bin) {
      splitbit val = flip ? ~s.state[i][bin] : s.state[i][bin];
      if (bin == last_bin) val &= last_mask;
      tmp[static_cast<size_t>(i) * n_bins + bin] = val;
    }
  }

  // --- 2. Sort split indices by canonical form ---
  const splitbit* t = tmp.data();
  std::iota(out.idx.begin(), out.idx.end(), split_int(0));
  std::sort(out.idx.begin(), out.idx.end(),
            [t, n_bins](split_int i, split_int j) {
              const splitbit* pi = t + static_cast<size_t>(i) * n_bins;
              const splitbit* pj = t + static_cast<size_t>(j) * n_bins;
              for (split_int bin = 0; bin < n_bins; ++bin) {
                if (pi[bin] < pj[bin]) return true;
                if (pi[bin] > pj[bin]) return false;
              }
              return false; // #nocov
            });

  // --- 3. Materialise canonical forms in sorted order ---
  for (split_int pos = 0; pos < n; ++pos) {
    std::memcpy(&out.canon[static_cast<size_t>(pos) * n_bins],
                t + static_cast<size_t>(out.idx[pos]) * n_bins,
                static_cast<size_t>(n_bins) * sizeof(splitbit));
  }
}

// Build canonical forms for a whole tree set.  Pure C++ (no R SEXP access), so
// safe to parallelise over trees.
static void build_canon_all(
    const std::vector<std::unique_ptr<SplitList>>& splits,
    const int32 n_tips,
    const int n_threads,
    std::vector<CanonSplits>& out
) {
  const int N = static_cast<int>(splits.size());
  out.resize(N);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_threads)
#else
  (void) n_threads;
#endif
  for (int k = 0; k < N; ++k) {
    build_canon(*splits[k], n_tips, out[k]);
  }
}

// ---------------------------------------------------------------------------
// O(n) exact-match detection for identical bipartitions.
//
// Merge-scans the two pre-sorted canonical split arrays built by build_canon().
//
// Returns the number of exact matches found.
// Writes results into scratch.a_match / scratch.b_match:
//   a_match[ai] = bi+1 if split ai matched split bi, else 0.
//   b_match[bi] = ai+1 if split bi matched split ai, else 0.
// ---------------------------------------------------------------------------
static split_int find_exact_matches(
    const CanonSplits& ca, const CanonSplits& cb,
    const split_int n_bins,
    MatchScratch& scratch
) {
  const split_int a_n = static_cast<split_int>(ca.idx.size());
  const split_int b_n = static_cast<split_int>(cb.idx.size());

  // Ensure buffers are large enough (grow lazily, never shrink)
  if (scratch.a_match.size() < static_cast<size_t>(a_n)) scratch.a_match.resize(a_n);
  if (scratch.b_match.size() < static_cast<size_t>(b_n)) scratch.b_match.resize(b_n);

  split_int* a_match = scratch.a_match.data();
  split_int* b_match = scratch.b_match.data();
  std::fill(a_match, a_match + a_n, split_int(0));
  std::fill(b_match, b_match + b_n, split_int(0));

  if (a_n == 0 || b_n == 0) return 0;

  const splitbit* ac = ca.canon.data();
  const splitbit* bc = cb.canon.data();
  const split_int* ai_of = ca.idx.data();
  const split_int* bi_of = cb.idx.data();

  split_int exact_n = 0;
  split_int ai_pos = 0, bi_pos = 0;
  while (ai_pos < a_n && bi_pos < b_n) {
    const splitbit* pa = ac + static_cast<size_t>(ai_pos) * n_bins;
    const splitbit* pb = bc + static_cast<size_t>(bi_pos) * n_bins;

    int cmp = 0;
    for (split_int bin = 0; bin < n_bins; ++bin) {
      if (pa[bin] < pb[bin]) { cmp = -1; break; }
      if (pa[bin] > pb[bin]) { cmp =  1; break; }
    }

    if (cmp < 0) {
      ++ai_pos;
    } else if (cmp > 0) {
      ++bi_pos;
    } else {
      const split_int ai = ai_of[ai_pos];
      const split_int bi = bi_of[bi_pos];
      a_match[ai] = bi + 1;
      b_match[bi] = ai + 1;
      ++exact_n;
      ++ai_pos;
      ++bi_pos;
    }
  }

  return exact_n;
}


template <typename Pop>
static double mutual_clustering_score_impl(
    const SplitList& a, const SplitList& b, const int32 n_tips,
    const CanonSplits& ca, const CanonSplits& cb,
    LapScratch& scratch, MatchScratch& mscratch
) {
  if (a.n_splits == 0 || b.n_splits == 0 || n_tips == 0) return 0.0;

  const split_int most_splits = std::max(a.n_splits, b.n_splits);
  const double n_tips_rcp = 1.0 / static_cast<double>(n_tips);

  constexpr cost max_score  = BIG;
  constexpr double over_max = 1.0 / static_cast<double>(BIG);
  const double max_over_tips = static_cast<double>(BIG) * n_tips_rcp;
  const double lg2_n = lg2_lookup(n_tips);

  // --- Phase 1: O(n log n) exact-match detection ---
  const split_int exact_n = find_exact_matches(ca, cb, a.n_bins, mscratch);
  const split_int* a_match = mscratch.a_match.data();
  const split_int* b_match = mscratch.b_match.data();

  // Accumulate exact-match score
  double exact_score = 0.0;
  for (split_int ai = 0; ai < a.n_splits; ++ai) {
    if (a_match[ai]) {
      const split_int na = a.in_split[ai];
      const split_int nA = n_tips - na;
      exact_score += TreeDist::ic_matching(na, nA, n_tips);
    }
  }

  // Early exit when everything matched exactly
  if (exact_n == b.n_splits || exact_n == a.n_splits) {
    return exact_score * n_tips_rcp;
  }

  // --- Phase 2: fill cost matrix for unmatched splits only (O(k²)) ---
  const split_int lap_n = most_splits - exact_n;

  // Build index maps for unmatched splits
  std::vector<split_int> a_unmatch, b_unmatch;
  a_unmatch.reserve(lap_n);
  b_unmatch.reserve(lap_n);
  for (split_int ai = 0; ai < a.n_splits; ++ai) {
    if (!a_match[ai]) a_unmatch.push_back(ai);
  }
  for (split_int bi = 0; bi < b.n_splits; ++bi) {
    if (!b_match[bi]) b_unmatch.push_back(bi);
  }

  scratch.score_pool.resize(lap_n);
  cost_matrix& score = scratch.score_pool;

  const split_int a_unmatched_n = static_cast<split_int>(a_unmatch.size());
  const split_int b_unmatched_n = static_cast<split_int>(b_unmatch.size());

  for (split_int a_pos = 0; a_pos < a_unmatched_n; ++a_pos) {
    const split_int ai   = a_unmatch[a_pos];
    const split_int na   = a.in_split[ai];
    const split_int nA   = n_tips - na;
    const auto* a_row = a.state[ai];

    const double offset_a = lg2_n - lg2_lookup(na);
    const double offset_A = lg2_n - lg2_lookup(nA);

    for (split_int b_pos = 0; b_pos < b_unmatched_n; ++b_pos) {
      const split_int bi   = b_unmatch[b_pos];
      const auto* b_row = b.state[bi];
      split_int a_and_b = 0;
      for (split_int bin = 0; bin < a.n_bins; ++bin) {
        a_and_b += Pop::count(a_row[bin] & b_row[bin]);
      }

      const split_int nb    = b.in_split[bi];
      const split_int nB    = n_tips - nb;
      const split_int a_and_B = na - a_and_b;
      const split_int A_and_b = nb - a_and_b;
      const split_int A_and_B = nA - A_and_b;

      if (a_and_b == A_and_b && a_and_b == a_and_B && a_and_b == A_and_B) {
        score(a_pos, b_pos) = max_score;
      } else {
        const double lg2_nb = lg2_lookup(nb);
        const double lg2_nB = lg2_lookup(nB);
        const double ic_sum =
          a_and_b * (lg2_lookup(a_and_b) + offset_a - lg2_nb) +
          a_and_B * (lg2_lookup(a_and_B) + offset_a - lg2_nB) +
          A_and_b * (lg2_lookup(A_and_b) + offset_A - lg2_nb) +
          A_and_B * (lg2_lookup(A_and_B) + offset_A - lg2_nB);
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

// Score-only version of mutual_clustering(). Thread-safe: uses only local
// storage and read-only globals (lg2 table, lookup tables).
// Passes allow_interrupt = false to lap() so it is safe to call from an
// OpenMP parallel region.
static double mutual_clustering_score(
    const SplitList& a, const SplitList& b, const int32 n_tips,
    const CanonSplits& ca, const CanonSplits& cb,
    LapScratch& scratch, MatchScratch& mscratch
) {
  return TreeDist::dispatch_popcount([&](auto pop) {
    return mutual_clustering_score_impl<decltype(pop)>(a, b, n_tips, ca, cb, scratch, mscratch);
  });
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
  TreeDist::check_ntip(n_tip);
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

  // Canonical forms depend on a single tree, so build them once per tree here
  // rather than once per pair inside find_exact_matches().
  std::vector<CanonSplits> canon;
  build_canon_all(splits, n_tip, n_threads, canon);

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
      res[p] = mutual_clustering_score(*splits[col], *splits[row], n_tip,
                                       canon[col], canon[row], scratch, mscratch);
    }
  }

  return result;
}


// =============================================================================
// InfoRobinsonFoulds (rf_info) — no LAP, thread-safe
// =============================================================================

template <typename Pop>
static double rf_info_score_impl(
    const SplitList& a, const SplitList& b, const int32 n_tips,
    const CanonSplits& ca, const CanonSplits& cb,
    MatchScratch& mscratch
) {
  const split_int a_n = a.n_splits;
  const split_int b_n = b.n_splits;
  if (a_n == 0 || b_n == 0) return 0;

  // Use sort+merge to find exact matches in O(n log n)
  const split_int exact_n = find_exact_matches(ca, cb, a.n_bins, mscratch);
  if (exact_n == 0) return 0;

  // Sum info contribution for each matched split in a
  const split_int* a_match = mscratch.a_match.data();
  const double lg2_unrooted_n = lg2_unrooted_lookup(n_tips);
  double score = 0;
  for (split_int ai = 0; ai < a_n; ++ai) {
    if (a_match[ai] == 0) continue;
    split_int leaves_in_split = 0;
    for (split_int bin = 0; bin < a.n_bins; ++bin) {
      leaves_in_split += Pop::count(a.state[ai][bin]);
    }
    score += lg2_unrooted_n
           - lg2_rooted_lookup(leaves_in_split)
           - lg2_rooted_lookup(n_tips - leaves_in_split);
  }
  return score;
}

static double rf_info_score(
    const SplitList& a, const SplitList& b, const int32 n_tips,
    const CanonSplits& ca, const CanonSplits& cb,
    MatchScratch& mscratch
) {
  return TreeDist::dispatch_popcount([&](auto pop) {
    return rf_info_score_impl<decltype(pop)>(a, b, n_tips, ca, cb, mscratch);
  });
}

//' @keywords internal
// [[Rcpp::export]]
NumericVector cpp_rf_info_all_pairs(
    const List& splits_list,
    const int   n_tip,
    const int   n_threads = 1
) {
  TreeDist::check_ntip(n_tip);
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

  std::vector<CanonSplits> canon;
  build_canon_all(splits, n_tip, n_threads, canon);

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
      res[p] = rf_info_score(*splits[col], *splits[row], n_tip,
                             canon[col], canon[row], mscratch);
    }
  }
  return result;
}


// =============================================================================
// MatchingSplitDistance — LAP-based
// =============================================================================

template <typename Pop>
static double msd_score_impl(
    const SplitList& a, const SplitList& b, const int32 n_tips,
    const CanonSplits& ca, const CanonSplits& cb,
    LapScratch& scratch, MatchScratch& mscratch
) {
  const split_int most_splits = std::max(a.n_splits, b.n_splits);
  if (most_splits == 0) return 0.0;
  const bool      a_has_more  = (a.n_splits > b.n_splits);
  const split_int a_extra     = a_has_more ? most_splits - b.n_splits : 0;
  const split_int b_extra     = a_has_more ? 0 : most_splits - a.n_splits;
  const split_int half_tips   = n_tips / 2;
  const cost      max_score   = BIG / most_splits;

  // --- Phase 1: O(n log n) exact-match detection ---
  const split_int exact_n = find_exact_matches(ca, cb, a.n_bins, mscratch);
  const split_int* a_match = mscratch.a_match.data();
  const split_int* b_match = mscratch.b_match.data();

  if (exact_n == b.n_splits || exact_n == a.n_splits) {
    return 0.0;
  }

  // --- Phase 2: fill cost matrix for unmatched splits only ---
  const split_int lap_n = most_splits - exact_n;

  std::vector<split_int> a_unmatch, b_unmatch;
  a_unmatch.reserve(lap_n);
  b_unmatch.reserve(lap_n);
  for (split_int ai = 0; ai < a.n_splits; ++ai) {
    if (!a_match[ai]) a_unmatch.push_back(ai);
  }
  for (split_int bi = 0; bi < b.n_splits; ++bi) {
    if (!b_match[bi]) b_unmatch.push_back(bi);
  }

  scratch.score_pool.resize(lap_n);
  cost_matrix& score = scratch.score_pool;

  const split_int a_unmatched_n = static_cast<split_int>(a_unmatch.size());
  const split_int b_unmatched_n = static_cast<split_int>(b_unmatch.size());

  for (split_int a_pos = 0; a_pos < a_unmatched_n; ++a_pos) {
    const split_int ai = a_unmatch[a_pos];
    for (split_int b_pos = 0; b_pos < b_unmatched_n; ++b_pos) {
      const split_int bi = b_unmatch[b_pos];
      splitbit total = 0;
      for (split_int bin = 0; bin < a.n_bins; ++bin) {
        total += Pop::count(a.state[ai][bin] ^ b.state[bi][bin]);
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

static double msd_score(
    const SplitList& a, const SplitList& b, const int32 n_tips,
    const CanonSplits& ca, const CanonSplits& cb,
    LapScratch& scratch, MatchScratch& mscratch
) {
  return TreeDist::dispatch_popcount([&](auto pop) {
    return msd_score_impl<decltype(pop)>(a, b, n_tips, ca, cb, scratch, mscratch);
  });
}

//' @keywords internal
// [[Rcpp::export]]
NumericVector cpp_msd_all_pairs(
    const List& splits_list,
    const int   n_tip,
    const int   n_threads = 1
) {
  TreeDist::check_ntip(n_tip);
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

  std::vector<CanonSplits> canon;
  build_canon_all(splits, n_tip, n_threads, canon);

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
      res[p] = msd_score(*splits[col], *splits[row], n_tip,
                         canon[col], canon[row], scratch, mscratch);
    }
  }
  return result;
}


// =============================================================================
// MatchingSplitInfo — LAP-based
// =============================================================================

template <typename Pop>
static double msi_score_impl(
    const SplitList& a, const SplitList& b, const int32 n_tips,
    LapScratch& scratch
) {
  const split_int most_splits = std::max(a.n_splits, b.n_splits);
  if (most_splits == 0) return 0.0;

  constexpr cost max_score = BIG;
  const double max_possible = lg2_unrooted_lookup(n_tips)
    - lg2_rooted_lookup((n_tips + 1) / 2)
    - lg2_rooted_lookup(n_tips / 2);
  const double score_over_possible = static_cast<double>(max_score) / max_possible;
  const double possible_over_score = max_possible / static_cast<double>(max_score);

  scratch.score_pool.resize(most_splits);
  cost_matrix& score = scratch.score_pool;

  for (split_int ai = 0; ai < a.n_splits; ++ai) {
    for (split_int bi = 0; bi < b.n_splits; ++bi) {
      split_int n_a_only = 0, n_a_and_b = 0, n_different = 0;
      splitbit different;
      for (split_int bin = 0; bin < a.n_bins; ++bin) {
        different   = a.state[ai][bin] ^ b.state[bi][bin];
        n_different += Pop::count(different);
        n_a_only   += Pop::count(a.state[ai][bin] &  different);
        n_a_and_b  += Pop::count(a.state[ai][bin] & ~different);
      }
      const split_int n_same = n_tips - n_different;
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

// NOTE: msi_score does NOT use exact-match detection.  For the MSI metric,
// matching a split to its identical copy is NOT necessarily globally optimal —
// the LAP can find a better assignment where the identical split pairs with
// a different (compatible, containing) split that has higher mmsi_score.
// This differs from MCI and MSD where exact matches are provably optimal.
static double msi_score(
    const SplitList& a, const SplitList& b, const int32 n_tips,
    LapScratch& scratch
) {
  return TreeDist::dispatch_popcount([&](auto pop) {
    return msi_score_impl<decltype(pop)>(a, b, n_tips, scratch);
  });
}

//' @keywords internal
// [[Rcpp::export]]
NumericVector cpp_msi_all_pairs(
    const List& splits_list,
    const int   n_tip,
    const int   n_threads = 1
) {
  TreeDist::check_ntip(n_tip);
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

template <bool Fast, typename Pop>
static inline void shared_phylo_fill_lap_impl(
    const SplitList& a, const SplitList& b, const int32 n_tips,
    const double best_overlap, const double score_over_possible,
    const cost max_score, cost_matrix& score
) {
  for (split_int ai = 0; ai < a.n_splits; ++ai) {
    for (split_int bi = 0; bi < b.n_splits; ++bi) {
      const double spi = TreeDist::spi_overlap<Fast, Pop>(
        a.state[ai], b.state[bi], n_tips,
        a.in_split[ai], b.in_split[bi], a.n_bins);
      score(ai, bi) = (spi < 0.0) ? max_score
                                  : cost((spi - best_overlap) * score_over_possible);
    }
    score.padRowAfterCol(ai, b.n_splits, max_score);
  }
}

// NOTE: For SPI, matching each split to its identical copy is NOT always
// globally optimal: the LAP can do better by pairing an identical split with
// a different compatible split.  shared_phylo_pinned() therefore treats
// identical matches as a guess that must be certified, and falls back to the
// full LAP when it cannot be.
// LAP-fill body for shared_phylo_score, templated on Fast so the per-cell
// spi_overlap call resolves at compile time to direct lg2_rooted[] access
// (Fast=true) or the bounds-checked lookup (Fast=false).
template<bool Fast>
static inline void shared_phylo_fill_lap(
    const SplitList& a, const SplitList& b, const int32 n_tips,
    const double best_overlap, const double score_over_possible,
    const cost max_score, cost_matrix& score
) {
  return TreeDist::dispatch_popcount([&](auto pop) {
    return shared_phylo_fill_lap_impl<Fast, decltype(pop)>(a, b, n_tips, best_overlap, score_over_possible, max_score, score);
  });
}

// Scaling shared by the full and pinned SPI solvers.
struct SpiScale {
  double best_overlap;
  double score_over_possible;
  double possible_over_score;
};

static SpiScale spi_scale(const int32 n_tips, const bool fast) {
  const split_int overlap_a = (n_tips + 1) / 2;
  double best_overlap, max_possible;
  if (fast) {
    best_overlap = TreeDist::one_overlap<true>(overlap_a, n_tips / 2, n_tips);
    max_possible = lg2_unrooted[n_tips] - best_overlap;
  } else {
    // LCOV_EXCL_START # nocov start
    best_overlap = TreeDist::one_overlap<false>(overlap_a, n_tips / 2, n_tips);
    max_possible = lg2_unrooted_lookup(n_tips) - best_overlap;
    // LCOV_EXCL_STOP # nocov end
  }
  return {best_overlap,
          static_cast<double>(BIG) / max_possible,
          max_possible / static_cast<double>(BIG)};
}

// Full SPI: fill the whole cost matrix and solve the LAP.
static double shared_phylo_full(
    const SplitList& a, const SplitList& b, const int32 n_tips,
    const bool fast, const SpiScale& sc, LapScratch& scratch
) {
  const split_int most_splits = std::max(a.n_splits, b.n_splits);
  constexpr cost max_score = BIG;

  scratch.score_pool.resize(most_splits);
  cost_matrix& score = scratch.score_pool;

  if (fast) {
    shared_phylo_fill_lap<true>(a, b, n_tips, sc.best_overlap,
                                sc.score_over_possible, max_score, score);
  } else {
    // LCOV_EXCL_START # nocov start
    shared_phylo_fill_lap<false>(a, b, n_tips, sc.best_overlap,
                                 sc.score_over_possible, max_score, score);
    // LCOV_EXCL_STOP # nocov end
  }
  score.padAfterRow(a.n_splits, max_score);

  scratch.ensure(most_splits);
  auto& rowsol = scratch.rowsol;
  auto& colsol = scratch.colsol;

  return static_cast<double>(
    (max_score * most_splits) - lap(most_splits, score, rowsol, colsol, false, scratch)
  ) * sc.possible_over_score;
}

template <typename Pop>
static double shared_phylo_pinned_impl(
    const SplitList& a, const SplitList& b, const int32 n_tips,
    const SpiScale& sc, const split_int* a_match,
    LapScratch& scratch, PinScratch& pin
) {
  const split_int most_splits = std::max(a.n_splits, b.n_splits);
  constexpr cost max_score = BIG;
  // Lower bound on the cost of matching two distinct compatible splits whose
  // smaller sides hold s and t tips: the cheaper of nesting the smaller in the
  // larger, or keeping them disjoint.  Less 1e-6 bits to absorb rounding.
  const double slack = 1e-6 * sc.score_over_possible + 2;
  auto bound_cost = [&](split_int s, split_int t) -> cost {
    if (s > t) std::swap(s, t);
    double ov = std::numeric_limits<double>::infinity();
    if (s < t) {
      ov = lg2_rooted[t] + lg2_rooted[n_tips - s] - lg2_rooted[t - s + 1];
    }
    if (s + t < n_tips) {
      ov = std::min(ov, lg2_rooted[n_tips - t] + lg2_rooted[n_tips - s] -
                          lg2_rooted[n_tips - t - s + 1]);
    }
    if (ov == std::numeric_limits<double>::infinity()) return max_score;
    return cost((ov - sc.best_overlap) * sc.score_over_possible - slack);
  };

  auto cell = [&](split_int ai, split_int bi) -> cost {
    const double spi = TreeDist::spi_overlap<true, Pop>(
      a.state[ai], b.state[bi], n_tips,
      a.in_split[ai], b.in_split[bi], a.n_bins);
    return (spi < 0.0) ? max_score
                       : cost((spi - sc.best_overlap) * sc.score_over_possible);
  };

  auto& pin_a  = pin.pin_a;
  auto& pin_b  = pin.pin_b;
  auto& pin_d  = pin.pin_d;
  auto& pin_lo = pin.pin_lo;
  auto& pin_u  = pin.pin_u;
  auto& free_a = pin.free_a;
  auto& free_b = pin.free_b;
  auto& u      = pin.u;
  pin_a.clear(); pin_b.clear(); pin_d.clear();
  free_a.clear(); free_b.clear();
  pin.b_pinned.assign(b.n_splits, 0);
  if (pin.group_of.size() < static_cast<size_t>(n_tips / 2 + 1)) {
    pin.group_of.assign(n_tips / 2 + 1, -1);
  }
  for (split_int ai = 0; ai < a.n_splits; ++ai) {
    if (a_match[ai]) {
      const split_int bi = a_match[ai] - 1;
      pin_a.push_back(ai);
      pin_b.push_back(bi);
      pin_d.push_back(cell(ai, bi));
      pin.b_pinned[bi] = 1;
    } else {
      free_a.push_back(ai);
    }
  }
  for (split_int bi = 0; bi < b.n_splits; ++bi) {
    if (!pin.b_pinned[bi]) free_b.push_back(bi);
  }

  // Each round either certifies or releases at least one pair, so the loop
  // ends; in the worst case every pair is released and the residual LAP is
  // the full problem.
  for (;;) {
    const split_int n_pin = static_cast<split_int>(pin_a.size());
    const split_int dim = most_splits - n_pin;
    const split_int fa = static_cast<split_int>(free_a.size());
    const split_int fb = static_cast<split_int>(free_b.size());

    // --- Residual LAP over unpinned splits, padded with max_score ---
    cost residual = 0;
    cost_matrix& score = scratch.score_pool;
    if (dim > 0) {
      score.resize(dim);
      for (split_int r = 0; r < fa; ++r) {
        for (split_int c = 0; c < fb; ++c) {
          score(r, c) = cell(free_a[r], free_b[c]);
        }
        score.padRowAfterCol(r, fb, max_score);
      }
      score.padAfterRow(fa, max_score);
      scratch.ensure(dim);
      residual = lap(dim, score, scratch.rowsol, scratch.colsol, false,
                     scratch);

      // Row duals; complementary slackness confirms the residual is optimal.
      const cost* v = scratch.v.data();
      u.resize(dim);
      for (split_int r = 0; r < dim; ++r) {
        const cost* row = score.row(r);
        cost min_r = row[0] - v[0];
        for (split_int c = 1; c < dim; ++c) {
          min_r = std::min(min_r, row[c] - v[c]);
        }
        u[r] = min_r;
        const split_int c_r = scratch.rowsol[r];
        if (row[c_r] - v[c_r] != min_r) return -1.0;
      }
    }

    // --- Dual window [lo, hi] for each pinned pair's u_i ---
    const cost* v = scratch.v.data();
    pin_lo.resize(n_pin);
    pin_u.resize(n_pin);
    split_int kept = 0;
    for (split_int p = 0; p < n_pin; ++p) {
      const split_int ai = pin_a[p];
      const split_int bi = pin_b[p];
      const cost d = pin_d[p];
      // u_i <= hi: residual columns in row ai
      cost hi = std::numeric_limits<cost>::max();
      for (split_int c = 0; c < fb; ++c) {
        hi = std::min(hi, cell(ai, free_b[c]) - v[c]);
      }
      for (split_int c = fb; c < dim; ++c) hi = std::min(hi, max_score - v[c]);
      // v_j = d - u_i <= residual rows in column bi
      cost lo = std::numeric_limits<cost>::min();
      for (split_int r = 0; r < fa; ++r) {
        lo = std::max(lo, d - (cell(free_a[r], bi) - u[r]));
      }
      for (split_int r = fa; r < dim; ++r) lo = std::max(lo, d - (max_score - u[r]));

      if (lo > hi) {
        free_a.push_back(ai);
        free_b.push_back(bi);
      } else {
        pin_a[kept] = ai;
        pin_b[kept] = bi;
        pin_d[kept] = d;
        pin_lo[kept] = lo;
        pin_u[kept] = hi;
        ++kept;
      }
    }
    if (kept < n_pin) {
      pin_a.resize(kept); pin_b.resize(kept); pin_d.resize(kept);
      continue;
    }

    // --- Cells between pinned pairs ---
    // With v_q = d_q - u_q, cell (p, q) needs u_p - w_q <= c(p, q), where
    // w_q = u_q - d_q.  Both splits come from tree a, so they are compatible,
    // and c(p, q) >= bound_cost(size_p, size_q).  Taking each u at the top of
    // its window (u = hi) keeps v small.  For each pair of size groups, walk
    // p by decreasing u and q by increasing w, evaluating cells only while
    // the bound fails to clear them.
    if (n_pin > 1) {
      auto& cls = pin.cls;
      auto& group_size = pin.group_size;
      auto& group_start = pin.group_start;
      group_size.clear();
      cls.resize(n_pin);
      for (split_int p = 0; p < n_pin; ++p) {
        const split_int in = a.in_split[pin_a[p]];
        const split_int size = std::min(in, split_int(n_tips - in));
        if (pin.group_of[size] < 0) {
          pin.group_of[size] = static_cast<split_int>(group_size.size());
          group_size.push_back(size);
        }
        cls[p] = pin.group_of[size];
      }
      for (const split_int size : group_size) pin.group_of[size] = -1;
      const split_int n_groups = static_cast<split_int>(group_size.size());
      group_start.assign(n_groups + 1, 0);
      for (split_int p = 0; p < n_pin; ++p) ++group_start[cls[p] + 1];
      for (split_int g = 0; g < n_groups; ++g) {
        group_start[g + 1] += group_start[g];
      }

      auto& by_u = pin.by_u;
      auto& by_w = pin.by_w;
      by_u.resize(n_pin);
      std::iota(by_u.begin(), by_u.end(), split_int(0));
      by_w = by_u;
      std::sort(by_u.begin(), by_u.end(), [&](split_int x, split_int y) {
        return cls[x] != cls[y] ? cls[x] < cls[y] : pin_u[x] > pin_u[y];
      });
      std::sort(by_w.begin(), by_w.end(), [&](split_int x, split_int y) {
        return cls[x] != cls[y] ? cls[x] < cls[y]
                                : pin_u[x] - pin_d[x] < pin_u[y] - pin_d[y];
      });

      for (split_int gp = 0; gp < n_groups; ++gp) {
        for (split_int gq = 0; gq < n_groups; ++gq) {
          const cost bound = bound_cost(group_size[gp], group_size[gq]);
          const split_int q_first = group_start[gq];
          const split_int q_end = group_start[gq + 1];
          const split_int q0 = by_w[q_first];
          for (split_int pi = group_start[gp]; pi < group_start[gp + 1];
               ++pi) {
            const split_int p = by_u[pi];
            if (pin_u[p] - (pin_u[q0] - pin_d[q0]) <= bound) break;
            for (split_int qi = q_first; qi < q_end; ++qi) {
              const split_int q = by_w[qi];
              const cost lhs = pin_u[p] - (pin_u[q] - pin_d[q]);
              if (lhs <= bound) break;
              if (p != q && lhs > cell(pin_a[p], pin_b[q])) {
                return -1.0; // # nocov: not observed; defensive
              }
            }
          }
        }
      }
    }

    cost total = residual;
    for (split_int p = 0; p < n_pin; ++p) total += pin_d[p];
    return static_cast<double>((max_score * most_splits) - total) *
      sc.possible_over_score;
  }
}

// Pinned SPI: match identical splits to one another, solve the LAP on the
// remaining splits only, then prove the combined assignment optimal.
//
// Matching identical splits is not always optimal (see NOTE above), so the
// pinned assignment is certified with LP duals (u_i + v_j <= c_ij for every
// cell, with equality on assigned cells).  The residual LAP supplies duals
// for unpinned rows and columns; each pinned pair (i, j) with cost d needs
// u_i + v_j = d, and its window for u_i is set by the residual cells in row
// i and column j.  A pair whose window is empty is released into the
// residual problem and the round repeats.
//
// Cells between two pinned pairs are never computed.  Pinned splits all
// belong to tree a, so any two are compatible, and the cost of matching them
// is bounded below by a function of their sizes alone; checking duals
// against that bound, one group of equal-sized splits at a time, covers
// every such cell.
//
// Returns a negative value if no certificate is found; the caller then solves
// the full problem.
static double shared_phylo_pinned(
    const SplitList& a, const SplitList& b, const int32 n_tips,
    const SpiScale& sc, const split_int* a_match,
    LapScratch& scratch, PinScratch& pin
) {
  return TreeDist::dispatch_popcount([&](auto pop) {
    return shared_phylo_pinned_impl<decltype(pop)>(a, b, n_tips, sc, a_match, scratch, pin);
  });
}

static double shared_phylo_score(
    const SplitList& a, const SplitList& b, const int32 n_tips,
    const CanonSplits& ca, const CanonSplits& cb,
    LapScratch& scratch, MatchScratch& mscratch, PinScratch& pin
) {
  const split_int most_splits = std::max(a.n_splits, b.n_splits);
  if (most_splits == 0) return 0.0;

  // Hot path: n_tips <= SL_MAX_TIPS, so all lg2_rooted/lg2_unrooted indices
  // fall in the precomputed table.  Dispatch once here so the per-cell loop
  // emits direct table loads with no bounds-check branch.
  // The slow branch is reachable only on TT 2.3 builds with n_tips > 32705,
  // which is too large to exercise in routine tests; the equivalent slow
  // template instantiations (one_overlap<false> etc.) are covered by the
  // legacy single-pair entry (cpp_shared_phylo) via force_slow=true.
  const bool fast = n_tips <= static_cast<int32>(SL_MAX_TIPS + 1);
  const SpiScale sc = spi_scale(n_tips, fast);

  if (fast) {
    const split_int exact_n = find_exact_matches(ca, cb, a.n_bins, mscratch);
    // Pinning pays only when most splits are shared.
    if (exact_n * 10 >= most_splits * 7) {
      const double pinned = shared_phylo_pinned(
        a, b, n_tips, sc, mscratch.a_match.data(), scratch, pin);
      if (pinned >= 0.0) return pinned;
    }
  }
  return shared_phylo_full(a, b, n_tips, fast, sc, scratch);
}

//' @keywords internal
// [[Rcpp::export]]
NumericVector cpp_shared_phylo_all_pairs(
    const List& splits_list,
    const int   n_tip,
    const int   n_threads = 1
) {
  TreeDist::check_ntip(n_tip);
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

  std::vector<CanonSplits> canon;
  build_canon_all(splits, n_tip, n_threads, canon);

  NumericVector result(n_pairs);
  double* const res = result.begin();

  const int n_scratch = std::max(1, n_threads);
  std::vector<LapScratch> scratches(n_scratch);
  std::vector<MatchScratch> mscratches(n_scratch);
  std::vector<PinScratch> pscratches(n_scratch);

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
    for (int row = col + 1; row < N; ++row) {
      const int p = col * (N - 1) - col * (col - 1) / 2 + row - col - 1;
      res[p] = shared_phylo_score(*splits[col], *splits[row], n_tip,
                                  canon[col], canon[row], scratches[tid],
                                  mscratches[tid], pscratches[tid]);
    }
  }
  return result;
}


// =============================================================================
// Jaccard / Nye similarity — LAP-based; k and allow_conflict are per-call
// =============================================================================

template <typename Pop>
static double jaccard_score_impl(
    const SplitList& a, const SplitList& b, const int32 n_tips,
    const CanonSplits& ca, const CanonSplits& cb,
    const double exponent, const bool allow_conflict,
    LapScratch& scratch, MatchScratch& mscratch
) {
  const split_int most_splits = std::max(a.n_splits, b.n_splits);
  if (most_splits == 0) return 0.0;

  constexpr cost   max_score  = BIG;
  constexpr double max_scoreL = static_cast<double>(max_score);

  // --- Phase 1: O(n log n) exact-match detection ---
  // Only used when allow_conflict=true; otherwise the full LAP may reassign
  // non-matching splits to compatible (non-exact) partners.
  split_int exact_n = 0;
  if (allow_conflict) {
    exact_n = find_exact_matches(ca, cb, a.n_bins, mscratch);
  } else {
    // Ensure match arrays are sized and zeroed
    if (mscratch.a_match.size() < static_cast<size_t>(a.n_splits))
      mscratch.a_match.resize(a.n_splits);
    if (mscratch.b_match.size() < static_cast<size_t>(b.n_splits))
      mscratch.b_match.resize(b.n_splits);
    std::fill(mscratch.a_match.data(), mscratch.a_match.data() + a.n_splits, split_int(0));
    std::fill(mscratch.b_match.data(), mscratch.b_match.data() + b.n_splits, split_int(0));
  }
  const split_int* a_match = mscratch.a_match.data();
  const split_int* b_match = mscratch.b_match.data();

  if (exact_n == b.n_splits || exact_n == a.n_splits) {
    return static_cast<double>(exact_n);
  }

  // --- Phase 2: fill cost matrix for unmatched splits only ---
  const split_int lap_n = most_splits - exact_n;

  std::vector<split_int> a_unmatch, b_unmatch;
  a_unmatch.reserve(lap_n);
  b_unmatch.reserve(lap_n);
  for (split_int ai = 0; ai < a.n_splits; ++ai) {
    if (!a_match[ai]) a_unmatch.push_back(ai);
  }
  for (split_int bi = 0; bi < b.n_splits; ++bi) {
    if (!b_match[bi]) b_unmatch.push_back(bi);
  }

  scratch.score_pool.resize(lap_n);
  cost_matrix& score = scratch.score_pool;

  const split_int a_unmatched_n = static_cast<split_int>(a_unmatch.size());
  const split_int b_unmatched_n = static_cast<split_int>(b_unmatch.size());

  for (split_int a_pos = 0; a_pos < a_unmatched_n; ++a_pos) {
    const split_int ai = a_unmatch[a_pos];
    const split_int na = a.in_split[ai];
    const split_int nA = n_tips - na;

    for (split_int b_pos = 0; b_pos < b_unmatched_n; ++b_pos) {
      const split_int bi = b_unmatch[b_pos];
      split_int a_and_b = 0;
      for (split_int bin = 0; bin < a.n_bins; ++bin) {
        a_and_b += Pop::count(a.state[ai][bin] & b.state[bi][bin]);
      }
      const split_int nb    = b.in_split[bi];
      const split_int nB    = n_tips - nb;
      const split_int a_and_B = na - a_and_b;
      const split_int A_and_b = nb - a_and_b;
      const split_int A_and_B = nB - a_and_B;

      if (!allow_conflict && !(
            a_and_b == na || a_and_B == na ||
            A_and_b == nA || A_and_B == nA)) {
        score(a_pos, b_pos) = max_score;
      } else {
        const split_int A_or_b = n_tips - a_and_B;
        const split_int a_or_B = n_tips - A_and_b;
        const split_int a_or_b = n_tips - A_and_B;
        const split_int A_or_B = n_tips - a_and_b;
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

static double jaccard_score(
    const SplitList& a, const SplitList& b, const int32 n_tips,
    const CanonSplits& ca, const CanonSplits& cb,
    const double exponent, const bool allow_conflict,
    LapScratch& scratch, MatchScratch& mscratch
) {
  return TreeDist::dispatch_popcount([&](auto pop) {
    return jaccard_score_impl<decltype(pop)>(a, b, n_tips, ca, cb, exponent, allow_conflict, scratch, mscratch);
  });
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
  TreeDist::check_ntip(n_tip);
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

  std::vector<CanonSplits> canon;
  build_canon_all(splits, n_tip, n_threads, canon);

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
      res[p] = jaccard_score(*splits[col], *splits[row], n_tip,
                             canon[col], canon[row], k, allow_conflict,
                             scratch, mscratch);
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
  TreeDist::check_ntip(n_tip);
  const int nA = splits_a.size();
  const int nB = splits_b.size();
  if (nA == 0 || nB == 0) return NumericMatrix(nA, nB);

  std::vector<std::unique_ptr<SplitList>> sa, sb;
  parse_split_list(splits_a, sa);
  parse_split_list(splits_b, sb);

  std::vector<CanonSplits> canon_a, canon_b;
  build_canon_all(sa, n_tip, n_threads, canon_a);
  build_canon_all(sb, n_tip, n_threads, canon_b);

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
    res[idx] = mutual_clustering_score(*sa[i], *sb[j], n_tip,
                                       canon_a[i], canon_b[j], scratch, mscratch);
  }
  return result;
}

//' @keywords internal
// [[Rcpp::export]]
NumericMatrix cpp_rf_info_cross_pairs(
    const List& splits_a, const List& splits_b,
    const int n_tip, const int n_threads = 1
) {
  TreeDist::check_ntip(n_tip);
  const int nA = splits_a.size();
  const int nB = splits_b.size();
  if (nA == 0 || nB == 0) return NumericMatrix(nA, nB);

  std::vector<std::unique_ptr<SplitList>> sa, sb;
  parse_split_list(splits_a, sa);
  parse_split_list(splits_b, sb);

  std::vector<CanonSplits> canon_a, canon_b;
  build_canon_all(sa, n_tip, n_threads, canon_a);
  build_canon_all(sb, n_tip, n_threads, canon_b);

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
    res[idx] = rf_info_score(*sa[i], *sb[j], n_tip,
                             canon_a[i], canon_b[j], mscratch);
  }
  return result;
}

//' @keywords internal
// [[Rcpp::export]]
NumericMatrix cpp_msd_cross_pairs(
    const List& splits_a, const List& splits_b,
    const int n_tip, const int n_threads = 1
) {
  TreeDist::check_ntip(n_tip);
  const int nA = splits_a.size();
  const int nB = splits_b.size();
  if (nA == 0 || nB == 0) return NumericMatrix(nA, nB);

  std::vector<std::unique_ptr<SplitList>> sa, sb;
  parse_split_list(splits_a, sa);
  parse_split_list(splits_b, sb);

  std::vector<CanonSplits> canon_a, canon_b;
  build_canon_all(sa, n_tip, n_threads, canon_a);
  build_canon_all(sb, n_tip, n_threads, canon_b);

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
    res[idx] = msd_score(*sa[i], *sb[j], n_tip,
                         canon_a[i], canon_b[j], scratch, mscratch);
  }
  return result;
}

//' @keywords internal
// [[Rcpp::export]]
NumericMatrix cpp_msi_cross_pairs(
    const List& splits_a, const List& splits_b,
    const int n_tip, const int n_threads = 1
) {
  TreeDist::check_ntip(n_tip);
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
  TreeDist::check_ntip(n_tip);
  const int nA = splits_a.size();
  const int nB = splits_b.size();
  if (nA == 0 || nB == 0) return NumericMatrix(nA, nB);

  std::vector<std::unique_ptr<SplitList>> sa, sb;
  parse_split_list(splits_a, sa);
  parse_split_list(splits_b, sb);

  std::vector<CanonSplits> canon_a, canon_b;
  build_canon_all(sa, n_tip, n_threads, canon_a);
  build_canon_all(sb, n_tip, n_threads, canon_b);

  NumericMatrix result(nA, nB);
  double* res = result.begin();

  const int n_scratch = std::max(1, n_threads);
  std::vector<LapScratch> scratches(n_scratch);
  std::vector<MatchScratch> mscratches(n_scratch);
  std::vector<PinScratch> pscratches(n_scratch);
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
    const int i = idx % nA;
    const int j = idx / nA;
    res[idx] = shared_phylo_score(*sa[i], *sb[j], n_tip, canon_a[i], canon_b[j],
                                  scratches[tid], mscratches[tid],
                                  pscratches[tid]);
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
  TreeDist::check_ntip(n_tip);
  const int nA = splits_a.size();
  const int nB = splits_b.size();
  if (nA == 0 || nB == 0) return NumericMatrix(nA, nB);

  std::vector<std::unique_ptr<SplitList>> sa, sb;
  parse_split_list(splits_a, sa);
  parse_split_list(splits_b, sb);

  std::vector<CanonSplits> canon_a, canon_b;
  build_canon_all(sa, n_tip, n_threads, canon_a);
  build_canon_all(sb, n_tip, n_threads, canon_b);

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
    res[idx] = jaccard_score(*sa[i], *sb[j], n_tip,
                             canon_a[i], canon_b[j], k, allow_conflict,
                             scratch, mscratch);
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
  TreeDist::check_ntip(n_tip);
  const int N = splits_list.size();
  NumericVector result(N);
  if (N == 0 || n_tip <= 0) return result;

  const double invN = 1.0 / static_cast<double>(n_tip);
  constexpr double invLog2 = 1.442695040888963387005;

  for (int i = 0; i < N; ++i) {
    SplitList sl(Rcpp::as<RawMatrix>(splits_list[i]));
    double total = 0.0;
    for (split_int s = 0; s < sl.n_splits; ++s) {
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
  TreeDist::check_ntip(n_tip);
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
    for (split_int s = 0; s < sl.n_splits; ++s) {
      const int k = sl.in_split[s];
      if (k < 2 || (n_tip - k) < 2) continue;
      total += l2u_n - l2r[k] - l2r[n_tip - k];
    }
    result[i] = total;
  }
  return result;
}


// =============================================================================
// Pairs of trees with different leaf sets
// =============================================================================
//
// Each pair is compared on the leaves the two trees share.  Restricting a
// tree's splits to a leaf subset, then discarding trivial splits and merging
// splits that become identical, gives the splits of the tree pruned to that
// subset; this avoids materialising pruned trees in R, one pair at a time.
//
// Splits arrive on each tree's own leaf order, with `tips[k]` giving the
// position of each of tree k's leaves in a label set shared by all trees.
// Shared leaves are numbered in the order they occur in the first tree of the
// pair, so both pruned trees index the same leaf with the same bit.

enum PrunedMetric {
  PM_NONE = 0, PM_MCI = 1, PM_IRF = 2, PM_MSD = 3, PM_MSI = 4, PM_SPI = 5,
  PM_JACCARD = 6, PM_RF = 7
};
enum PrunedInfo {
  PI_NONE = 0, PI_CLUSTERING = 1, PI_SPLITWISE = 2, PI_NSPLITS = 3
};

// A tree pruned to the leaves it shares with its partner.
//
// SplitList can only be built from an R RawMatrix, which is not thread-safe, so
// each thread builds one up front, sized for the largest input, and overwrites
// its public fields for every pair.  This relies on the row pointers in
// `state` having a stride of at least the pruned number of bins, which holds
// because pruning never adds splits or leaves.
struct PrunedTree {
  std::unique_ptr<SplitList> sl;
  CanonSplits canon;
  std::vector<splitbit> buf;     // pruned splits before deduplication
  std::vector<split_int> order;
};

struct PrunedScratch {
  PrunedTree a, b;
  std::vector<int32> rank;       // indexed by label: shared-leaf number, or < 0
  std::vector<int32> local_a, local_b;  // shared-leaf number of each leaf, or -1
  LapScratch lap;
  MatchScratch match;
  PinScratch pin;
};

// Restrict `src` to the shared leaves numbered by `local` (-1 if not shared),
// dropping trivial and duplicate splits; write the result to `out`.
static void prune_tree(const SplitList& src, const std::vector<int32>& local,
                       const int32 m, PrunedTree& out) {
  SplitList& dst = *out.sl;
  const int32 n_bins = std::max(int32(1), (m + SL_BIN_SIZE - 1) / SL_BIN_SIZE);
  const int32 n_local = static_cast<int32>(local.size());
  const splitbit last_mask = (m % SL_BIN_SIZE == 0)
    ? ~splitbit(0)
    : (splitbit(1) << (m % SL_BIN_SIZE)) - 1;

  out.buf.assign(static_cast<size_t>(src.n_splits) * n_bins, 0);
  split_int kept = 0;
  for (split_int s = 0; s < src.n_splits; ++s) {
    splitbit* row = &out.buf[static_cast<size_t>(kept) * n_bins];
    const splitbit* in = src.state[s];
    int32 count = 0;
    for (int32 t = 0; t < n_local; ++t) {
      const int32 r = local[t];
      if (r >= 0 && ((in[t / SL_BIN_SIZE] >> (t % SL_BIN_SIZE)) & 1)) {
        row[r / SL_BIN_SIZE] |= splitbit(1) << (r % SL_BIN_SIZE);
        ++count;
      }
    }
    if (count < 2 || m - count < 2) {
      std::fill(row, row + n_bins, splitbit(0));
      continue;
    }
    // Canonical form: shared leaf 0 on the "in" side
    if (!(row[0] & 1)) {
      for (int32 bin = 0; bin < n_bins; ++bin) row[bin] = ~row[bin];
      row[n_bins - 1] &= last_mask;
    }
    ++kept;
  }

  const splitbit* t = out.buf.data();
  out.order.resize(kept);
  std::iota(out.order.begin(), out.order.end(), split_int(0));
  std::sort(out.order.begin(), out.order.end(),
            [t, n_bins](split_int i, split_int j) {
              const splitbit* pi = t + static_cast<size_t>(i) * n_bins;
              const splitbit* pj = t + static_cast<size_t>(j) * n_bins;
              for (int32 bin = 0; bin < n_bins; ++bin) {
                if (pi[bin] != pj[bin]) return pi[bin] < pj[bin];
              }
              return false;
            });

  // Sorted, so duplicates are adjacent
  std::vector<splitbit>& canon = out.canon.canon;
  canon.resize(static_cast<size_t>(kept) * n_bins);
  split_int n_unique = 0;
  for (split_int pos = 0; pos < kept; ++pos) {
    const splitbit* row = t + static_cast<size_t>(out.order[pos]) * n_bins;
    splitbit* next = canon.data() + static_cast<size_t>(n_unique) * n_bins;
    if (n_unique > 0 && std::equal(row, row + n_bins, next - n_bins)) continue;
    std::memcpy(next, row, static_cast<size_t>(n_bins) * sizeof(splitbit));
    ++n_unique;
  }
  canon.resize(static_cast<size_t>(n_unique) * n_bins);
  out.canon.idx.resize(n_unique);
  std::iota(out.canon.idx.begin(), out.canon.idx.end(), split_int(0));

  dst.n_splits = n_unique;
  dst.n_bins = n_bins;
  for (split_int s = 0; s < n_unique; ++s) {
    const splitbit* row = &canon[static_cast<size_t>(s) * n_bins];
    int32 count = 0;
    for (int32 bin = 0; bin < n_bins; ++bin) {
      dst.state[s][bin] = row[bin];
      count += TreeTools::count_bits(row[bin]);
    }
    dst.in_split[s] = count;
  }
}

// Information content of a pruned tree on `m` leaves; matches
// cpp_clustering_entropy_batch() / cpp_splitwise_info_batch() / NSplits().
static double pruned_info(const SplitList& s, const int32 m, const int info,
                          const std::vector<double>& l2r) {
  switch (info) {
  case PI_CLUSTERING: {
    constexpr double invLog2 = 1.442695040888963387005;
    const double invN = 1.0 / static_cast<double>(m);
    double total = 0.0;
    for (split_int i = 0; i < s.n_splits; ++i) {
      const double p = s.in_split[i] * invN;
      total += -(p * std::log(p) + (1.0 - p) * std::log1p(-p));
    }
    return total * invLog2;
  }
  case PI_SPLITWISE: {
    if (s.n_splits == 0) return 0.0;
    const double l2u_m = l2r[m] - std::log2(2 * m - 3);
    double total = 0.0;
    for (split_int i = 0; i < s.n_splits; ++i) {
      const int32 k = s.in_split[i];
      total += l2u_m - l2r[k] - l2r[m - k];
    }
    return total;
  }
  case PI_NSPLITS:
    return static_cast<double>(s.n_splits);
  default: return NA_REAL; // #nocov
  }
}

static void pruned_pair(const SplitList& sa, const std::vector<int32>& tips_a,
                        const SplitList& sb, const std::vector<int32>& tips_b,
                        const int metric, const int info,
                        const double k, const bool allow_conflict,
                        const std::vector<double>& l2r, PrunedScratch& sc,
                        double& score, double& info_a, double& info_b) {
  // Flag b's leaves, then number the shared leaves in a's order
  constexpr int32 IN_B = -2;
  for (const int32 u : tips_b) sc.rank[u] = IN_B;
  const int32 n_a = static_cast<int32>(tips_a.size());
  sc.local_a.resize(n_a);
  int32 m = 0;
  for (int32 t = 0; t < n_a; ++t) {
    int32& r = sc.rank[tips_a[t]];
    if (r == IN_B) {
      r = m++;
      sc.local_a[t] = r;
    } else {
      sc.local_a[t] = -1;
    }
  }
  const int32 n_b = static_cast<int32>(tips_b.size());
  sc.local_b.resize(n_b);
  for (int32 t = 0; t < n_b; ++t) {
    const int32 r = sc.rank[tips_b[t]];
    sc.local_b[t] = r >= 0 ? r : -1;
  }
  // Only b's leaves were touched
  for (const int32 u : tips_b) sc.rank[u] = -1;

  prune_tree(sa, sc.local_a, m, sc.a);
  prune_tree(sb, sc.local_b, m, sc.b);
  const SplitList& a = *sc.a.sl;
  const SplitList& b = *sc.b.sl;

  if (info != PI_NONE) {
    info_a = pruned_info(a, m, info, l2r);
    info_b = pruned_info(b, m, info, l2r);
  }

  if (metric == PM_NONE) return;
  if (m < 4) {
    score = 0.0;
    return;
  }
  const CanonSplits& ca = sc.a.canon;
  const CanonSplits& cb = sc.b.canon;
  switch (metric) {
  case PM_MCI:
    score = mutual_clustering_score(a, b, m, ca, cb, sc.lap, sc.match);
    break;
  case PM_IRF:
    score = rf_info_score(a, b, m, ca, cb, sc.match);
    break;
  case PM_MSD:
    score = msd_score(a, b, m, ca, cb, sc.lap, sc.match);
    break;
  case PM_MSI:
    score = msi_score(a, b, m, sc.lap);
    break;
  case PM_SPI:
    score = shared_phylo_score(a, b, m, ca, cb, sc.lap, sc.match, sc.pin);
    break;
  case PM_JACCARD:
    score = jaccard_score(a, b, m, ca, cb, k, allow_conflict, sc.lap, sc.match);
    break;
  case PM_RF:
    score = a.n_splits + b.n_splits -
      2.0 * find_exact_matches(ca, cb, a.n_bins, sc.match);
    break;
  default: score = NA_REAL; // #nocov
  }
}

static void parse_tips(const List& list, const int n_union,
                       std::vector<std::vector<int32>>& out) {
  const int n = list.size();
  out.resize(n);
  for (int i = 0; i < n; ++i) {
    const Rcpp::IntegerVector v(list[i]);
    out[i].assign(v.begin(), v.end());
    for (const int32 u : out[i]) {
      if (u < 0 || u >= n_union) {
        Rcpp::stop("Leaf index out of range");                    // #nocov
      }
    }
  }
}

// Compare trees with different leaves, each pair on the leaves it shares.
//
// splits_a, splits_b: lists of raw split matrices, each on its tree's own leaf
//   order.  splits_b and tips_b are ignored if all_pairs.
// tips_a, tips_b: for each tree, the zero-based index of each of its leaves in
//   a label set of size n_union.
// all_pairs: compare all pairs within splits_a, in combn() order; otherwise
//   each tree in a with each tree in b, in column-major nA x nB order.
// metric: a PrunedMetric; info: a PrunedInfo, the information content of each
//   pruned tree to report.
// Returns list(score, info_a, info_b).
// [[Rcpp::export]]
List cpp_pruned_pairs(
    const List& splits_a, const List& tips_a,
    const List& splits_b, const List& tips_b,
    const bool all_pairs, const int n_union,
    const int metric, const int info,
    const double k = 1.0, const bool allow_conflict = true,
    const int n_threads = 1
) {
  TreeDist::check_ntip(n_union);
  std::vector<std::unique_ptr<SplitList>> sa, sb;
  std::vector<std::vector<int32>> ta, tb;
  parse_split_list(splits_a, sa);
  parse_tips(tips_a, n_union, ta);
  if (!all_pairs) {
    parse_split_list(splits_b, sb);
    parse_tips(tips_b, n_union, tb);
  }
  if (ta.size() != sa.size() || tb.size() != sb.size()) {
    Rcpp::stop("Each tree needs splits and leaf indices");        // #nocov
  }

  const int nA = static_cast<int>(sa.size());
  const int nB = all_pairs ? nA : static_cast<int>(sb.size());
  const R_xlen_t total = all_pairs
    ? static_cast<R_xlen_t>(nA) * (nA - 1) / 2
    : static_cast<R_xlen_t>(nA) * nB;

  NumericVector score(total, NA_REAL), info_a(total, NA_REAL),
    info_b(total, NA_REAL);
  if (total == 0) {
    return List::create(Rcpp::Named("score") = score,
                        Rcpp::Named("info_a") = info_a,
                        Rcpp::Named("info_b") = info_b);
  }

  int32 max_splits = 1, max_tips = 1;
  for (const auto* set : {&sa, &sb}) {
    for (const auto& s : *set) max_splits = std::max(max_splits, s->n_splits);
  }
  for (const auto* set : {&ta, &tb}) {
    for (const auto& t : *set) {
      max_tips = std::max(max_tips, static_cast<int32>(t.size()));
    }
  }

  // Log2Rooted(m) = log2((2m - 3)!!), as in cpp_splitwise_info_batch()
  std::vector<double> l2r(n_union + 1, 0.0);
  for (int m = 3; m <= n_union; ++m) {
    l2r[m] = l2r[m - 1] + std::log2(2 * m - 3);
  }

  const int n_scratch = std::max(1, n_threads);
  std::vector<PrunedScratch> scratch(n_scratch);
  {
    const RawMatrix blank(max_splits, (max_tips + R_BIN_SIZE - 1) / R_BIN_SIZE);
    for (PrunedScratch& sc : scratch) {
      sc.a.sl = std::make_unique<SplitList>(blank);
      sc.b.sl = std::make_unique<SplitList>(blank);
      sc.rank.assign(n_union, -1);
    }
  }

  double* const res = score.begin();
  double* const res_a = info_a.begin();
  double* const res_b = info_b.begin();
  const int n_outer = all_pairs ? nA - 1 : nB;

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(n_threads)
#endif
  for (int j = 0; j < n_outer; ++j) {
#ifndef _OPENMP
    Rcpp::checkUserInterrupt();
#endif
#ifdef _OPENMP
    const int tid = omp_get_thread_num();
#else
    const int tid = 0;
#endif
    PrunedScratch& sc = scratch[tid];
    if (all_pairs) {
      // Pair (j, i), j < i, is entry j(N-1) - j(j-1)/2 + i - j - 1 of combn()
      const R_xlen_t base = static_cast<R_xlen_t>(j) * (nA - 1) -
        static_cast<R_xlen_t>(j) * (j - 1) / 2 - j - 1;
      for (int i = j + 1; i < nA; ++i) {
        const R_xlen_t p = base + i;
        pruned_pair(*sa[j], ta[j], *sa[i], ta[i], metric, info, k,
                    allow_conflict, l2r, sc, res[p], res_a[p], res_b[p]);
      }
    } else {
      for (int i = 0; i < nA; ++i) {
        const R_xlen_t p = static_cast<R_xlen_t>(j) * nA + i;
        pruned_pair(*sa[i], ta[i], *sb[j], tb[j], metric, info, k,
                    allow_conflict, l2r, sc, res[p], res_a[p], res_b[p]);
      }
    }
  }

  return List::create(Rcpp::Named("score") = score,
                      Rcpp::Named("info_a") = info_a,
                      Rcpp::Named("info_b") = info_b);
}
