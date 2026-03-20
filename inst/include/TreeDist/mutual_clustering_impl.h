#ifndef TREEDIST_MUTUAL_CLUSTERING_IMPL_H_
#define TREEDIST_MUTUAL_CLUSTERING_IMPL_H_

// Mutual Clustering Information — implementation.
//
// Guard with #define TREEDIST_MCI_IMPLEMENTATION before including.
// Include in exactly one translation unit per package (the same TU
// that defines TREEDIST_LAP_IMPLEMENTATION).
//
// Usage:
//   #define TREEDIST_LAP_IMPLEMENTATION
//   #define TREEDIST_MCI_IMPLEMENTATION
//   #include <TreeDist/lap_impl.h>
//   #include <TreeDist/mutual_clustering_impl.h>

#ifdef TREEDIST_MCI_IMPLEMENTATION

#include "mutual_clustering.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

namespace TreeDist {

// ---- Table definitions ----

double lg2[SL_MAX_TIPS + 1];
double lg2_double_factorial[SL_MAX_TIPS + SL_MAX_TIPS - 2];
double lg2_unrooted[SL_MAX_TIPS + 2];
double* lg2_rooted = &lg2_unrooted[0] + 1;

void init_lg2_tables(int max_tips) {
  lg2[0] = 0;
  const int lg2_limit = std::min(max_tips + 1,
                                 static_cast<int>(SL_MAX_TIPS + 1));
  for (int i = 1; i < lg2_limit; ++i) {
    lg2[i] = std::log2(static_cast<double>(i));
  }

  for (int i = 0; i < 3; ++i) {
    lg2_double_factorial[i] = 0;
    lg2_unrooted[i] = 0;
  }

  const int df_limit = std::min(2 * max_tips,
                                static_cast<int>(SL_MAX_TIPS + SL_MAX_TIPS - 2));
  for (int i = 2; i < df_limit; ++i) {
    lg2_double_factorial[i] =
      lg2_double_factorial[i - 2] + std::log2(static_cast<double>(i));
  }

  const int ur_limit = std::min(max_tips + 2,
                                static_cast<int>(SL_MAX_TIPS + 2));
  for (int i = 3; i < ur_limit; ++i) {
    lg2_unrooted[i] = lg2_double_factorial[i + i - 5];
  }
}

// ---- Sort+merge exact-match detection (internal) ----
//
// Canonicalise each split so bit 0 is always set (flip complement if not),
// sort by canonical form, merge-scan.  O(n log n).
//
// Writes a_match[ai] = bi+1 if matched, else 0.  Likewise b_match.
// Returns number of exact matches.

namespace detail {

static int16 find_exact_matches_raw(
    const splitbit* const* a_state, const int16* /*a_in*/, int16 a_n,
    const splitbit* const* b_state, const int16* /*b_in*/, int16 b_n,
    int16 n_bins, int32 n_tips,
    int16* a_match, int16* b_match)
{
  std::fill(a_match, a_match + a_n, int16(0));
  std::fill(b_match, b_match + b_n, int16(0));
  if (a_n == 0 || b_n == 0) return 0;

  const int16 last_bin = n_bins - 1;
  const splitbit last_mask = (n_tips % SL_BIN_SIZE == 0)
    ? ~splitbit(0)
    : (splitbit(1) << (n_tips % SL_BIN_SIZE)) - 1;

  // Flat buffers for canonical forms
  std::vector<splitbit> a_canon(static_cast<std::size_t>(a_n) * n_bins);
  std::vector<splitbit> b_canon(static_cast<std::size_t>(b_n) * n_bins);

  for (int16 i = 0; i < a_n; ++i) {
    const bool flip = !(a_state[i][0] & 1);
    for (int16 bin = 0; bin < n_bins; ++bin) {
      splitbit val = flip ? ~a_state[i][bin] : a_state[i][bin];
      if (bin == last_bin) val &= last_mask;
      a_canon[i * n_bins + bin] = val;
    }
  }
  for (int16 i = 0; i < b_n; ++i) {
    const bool flip = !(b_state[i][0] & 1);
    for (int16 bin = 0; bin < n_bins; ++bin) {
      splitbit val = flip ? ~b_state[i][bin] : b_state[i][bin];
      if (bin == last_bin) val &= last_mask;
      b_canon[i * n_bins + bin] = val;
    }
  }

  // Sort index arrays by canonical form
  auto canon_less = [&](const splitbit* canon, int16 n_b, int16 i, int16 j) {
    for (int16 bin = 0; bin < n_b; ++bin) {
      const splitbit vi = canon[i * n_b + bin];
      const splitbit vj = canon[j * n_b + bin];
      if (vi < vj) return true;
      if (vi > vj) return false;
    }
    return false;
  };

  std::vector<int16> a_order(a_n), b_order(b_n);
  std::iota(a_order.begin(), a_order.end(), int16(0));
  std::iota(b_order.begin(), b_order.end(), int16(0));

  std::sort(a_order.begin(), a_order.end(),
            [&](int16 i, int16 j) {
              return canon_less(a_canon.data(), n_bins, i, j);
            });
  std::sort(b_order.begin(), b_order.end(),
            [&](int16 i, int16 j) {
              return canon_less(b_canon.data(), n_bins, i, j);
            });

  // Merge-scan
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

} // namespace detail


// ---- MCI score implementation ----

double mutual_clustering_score(
    const splitbit* const* a_state, const int16* a_in, int16 a_n_splits,
    const splitbit* const* b_state, const int16* b_in, int16 b_n_splits,
    int16 n_bins, int32 n_tips,
    LapScratch& scratch)
{
  if (a_n_splits == 0 || b_n_splits == 0 || n_tips == 0) return 0.0;

  const int16 most_splits = std::max(a_n_splits, b_n_splits);
  const double n_tips_rcp = 1.0 / static_cast<double>(n_tips);

  constexpr cost max_score  = BIG;
  constexpr double over_max = 1.0 / static_cast<double>(BIG);
  const double max_over_tips = static_cast<double>(BIG) * n_tips_rcp;
  const double lg2_n = lg2[n_tips];

  // --- Phase 1: O(n log n) exact-match detection ---
  std::vector<int16> a_match_buf(a_n_splits);
  std::vector<int16> b_match_buf(b_n_splits);

  const int16 exact_n = detail::find_exact_matches_raw(
    a_state, a_in, a_n_splits,
    b_state, b_in, b_n_splits,
    n_bins, n_tips,
    a_match_buf.data(), b_match_buf.data());

  const int16* a_match = a_match_buf.data();
  const int16* b_match = b_match_buf.data();

  // Accumulate exact-match score
  double exact_score = 0.0;
  for (int16 ai = 0; ai < a_n_splits; ++ai) {
    if (a_match[ai]) {
      const int16 na = a_in[ai];
      const int16 nA = static_cast<int16>(n_tips - na);
      exact_score += ic_matching(na, nA, static_cast<int16>(n_tips));
    }
  }

  // Early exit when everything matched exactly
  if (exact_n == b_n_splits || exact_n == a_n_splits) {
    return exact_score * n_tips_rcp;
  }

  // --- Phase 2: fill cost matrix for unmatched splits only (O(k²)) ---
  const int16 lap_n = most_splits - exact_n;

  std::vector<int16> a_unmatch, b_unmatch;
  a_unmatch.reserve(lap_n);
  b_unmatch.reserve(lap_n);
  for (int16 ai = 0; ai < a_n_splits; ++ai) {
    if (!a_match[ai]) a_unmatch.push_back(ai);
  }
  for (int16 bi = 0; bi < b_n_splits; ++bi) {
    if (!b_match[bi]) b_unmatch.push_back(bi);
  }

  scratch.score_pool.resize(lap_n);
  CostMatrix& score = scratch.score_pool;

  const int16 a_unmatched_n = static_cast<int16>(a_unmatch.size());
  const int16 b_unmatched_n = static_cast<int16>(b_unmatch.size());

  for (int16 a_pos = 0; a_pos < a_unmatched_n; ++a_pos) {
    const int16 ai = a_unmatch[a_pos];
    const int16 na = a_in[ai];
    const int16 nA = static_cast<int16>(n_tips - na);
    const splitbit* a_row = a_state[ai];

    const double offset_a = lg2_n - lg2[na];
    const double offset_A = lg2_n - lg2[nA];

    for (int16 b_pos = 0; b_pos < b_unmatched_n; ++b_pos) {
      const int16 bi = b_unmatch[b_pos];
      const splitbit* b_row = b_state[bi];
      int16 a_and_b = 0;
      for (int16 bin = 0; bin < n_bins; ++bin) {
        a_and_b += TreeTools::count_bits(a_row[bin] & b_row[bin]);
      }

      const int16 nb = b_in[bi];
      const int16 nB = static_cast<int16>(n_tips - nb);
      const int16 a_and_B = na - a_and_b;
      const int16 A_and_b = nb - a_and_b;
      const int16 A_and_B = nA - A_and_b;

      if (a_and_b == A_and_b && a_and_b == a_and_B
          && a_and_b == A_and_B) {
        score(a_pos, b_pos) = max_score;
      } else {
        const double lg2_nb = lg2[nb];
        const double lg2_nB = lg2[nB];
        const double ic_sum =
          a_and_b * (lg2[a_and_b] + offset_a - lg2_nb) +
          a_and_B * (lg2[a_and_B] + offset_a - lg2_nB) +
          A_and_b * (lg2[A_and_b] + offset_A - lg2_nb) +
          A_and_B * (lg2[A_and_B] + offset_A - lg2_nB);
        score(a_pos, b_pos) =
          max_score - static_cast<cost>(ic_sum * max_over_tips);
      }
    }
    if (b_unmatched_n < lap_n) {
      score.padRowAfterCol(a_pos, b_unmatched_n, max_score);
    }
  }
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

} // namespace TreeDist

#endif // TREEDIST_MCI_IMPLEMENTATION
#endif // TREEDIST_MUTUAL_CLUSTERING_IMPL_H_
