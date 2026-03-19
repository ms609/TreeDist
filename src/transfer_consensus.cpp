/* transfer_consensus.cpp
 *
 * C++ implementation of the Transfer Consensus algorithm (Takazawa 2025).
 * Replaces the R-level .PoolSplits, .TransferDistMat, .ComputeTD,
 * .CompatMat, and .GreedyBest/.GreedyFirst helpers.
 *
 * Exported function: cpp_transfer_consensus()
 */

#include <TreeTools/SplitList.h>
#include <Rcpp/Lightest>
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <functional>
#include <numeric>
#include <unordered_map>
#include <vector>

using Rcpp::IntegerVector;
using Rcpp::List;
using Rcpp::LogicalVector;
using Rcpp::RawMatrix;
using TreeTools::count_bits;
using TreeTools::SplitList;


// ============================================================================
// Types
// ============================================================================

using int16 = int_fast16_t;
using int32 = int_fast32_t;

// ============================================================================
// Pooled split representation
// ============================================================================

struct PooledSplits {
  int n_splits;
  int n_bins;
  int n_tips;
  splitbit last_mask;

  // Flat array of canonical split data: n_splits * n_bins elements.
  // Split i occupies [i*n_bins .. (i+1)*n_bins).
  std::vector<splitbit> data;

  // Per-split metadata
  std::vector<int> count;       // how many trees contain this split
  std::vector<int> light_side;  // min(popcount, n_tips - popcount)

  // Tree membership: tree_members[t] = list of unique-split indices in tree t
  std::vector<std::vector<int>> tree_members;

  // Access helpers
  const splitbit* split(int i) const { return &data[i * n_bins]; }
  splitbit*       split(int i)       { return &data[i * n_bins]; }
};


// ============================================================================
// FNV-1a hash for canonical split arrays
// ============================================================================

struct SplitHash {
  int n_bins;
  explicit SplitHash(int nb) : n_bins(nb) {}
  SplitHash() : n_bins(0) {}

  std::size_t operator()(const splitbit* sp) const {
    std::size_t h = 14695981039346656037ULL;
    const unsigned char* p = reinterpret_cast<const unsigned char*>(sp);
    const int bytes = n_bins * static_cast<int>(sizeof(splitbit));
    for (int i = 0; i < bytes; ++i) {
      h ^= static_cast<std::size_t>(p[i]);
      h *= 1099511628211ULL;
    }
    return h;
  }
};

struct SplitEqual {
  int n_bins;
  explicit SplitEqual(int nb) : n_bins(nb) {}
  SplitEqual() : n_bins(0) {}

  bool operator()(const splitbit* a, const splitbit* b) const {
    for (int bin = 0; bin < n_bins; ++bin) {
      if (a[bin] != b[bin]) return false;
    }
    return true;
  }
};

// ============================================================================
// pool_splits: deduplicate and canonicalise all splits from all trees
// ============================================================================

static PooledSplits pool_splits(const List& splits_list, int n_tips) {
  const int n_tree = splits_list.size();

  // Determine n_bins from first tree
  const RawMatrix first_mat = Rcpp::as<RawMatrix>(splits_list[0]);
  const int n_bytes = first_mat.ncol();
  const int n_bins = (n_bytes + static_cast<int>(sizeof(splitbit)) - 1)
    / static_cast<int>(sizeof(splitbit));

  const splitbit last_mask = (n_tips % SL_BIN_SIZE == 0)
    ? ~splitbit(0)
    : (splitbit(1) << (n_tips % SL_BIN_SIZE)) - 1;

  // Hash map from canonical split → unique index
  SplitHash hasher(n_bins);
  SplitEqual eq(n_bins);
  std::unordered_map<const splitbit*, int, SplitHash, SplitEqual>
    split_map(64, hasher, eq);

  // Temp buffer for canonicalisation
  std::vector<splitbit> canon_buf(n_bins);

  PooledSplits pool;
  pool.n_tips = n_tips;
  pool.n_bins = n_bins;
  pool.last_mask = last_mask;
  pool.n_splits = 0;
  pool.tree_members.resize(n_tree);

  for (int t = 0; t < n_tree; ++t) {
    SplitList sl(Rcpp::as<RawMatrix>(splits_list[t]));
    std::vector<int>& members = pool.tree_members[t];
    members.reserve(sl.n_splits);

    for (int16 s = 0; s < sl.n_splits; ++s) {
      // Copy split bits into canon_buf
      for (int bin = 0; bin < n_bins; ++bin) {
        canon_buf[bin] = sl.state[s][bin];
      }
      // Mask last bin
      canon_buf[n_bins - 1] &= last_mask;

      // Canonicalise: if bit 0 is set, flip
      if (canon_buf[0] & splitbit(1)) {
        for (int bin = 0; bin < n_bins; ++bin) {
          canon_buf[bin] = ~canon_buf[bin];
        }
        canon_buf[n_bins - 1] &= last_mask;
      }

      // Look up in map
      auto it = split_map.find(canon_buf.data());
      int idx;
      if (it != split_map.end()) {
        idx = it->second;
        pool.count[idx]++;
      } else {
        idx = pool.n_splits++;
        // Append canonical data
        const size_t old_sz = pool.data.size();
        pool.data.resize(old_sz + n_bins);
        std::copy(canon_buf.begin(), canon_buf.end(),
                  pool.data.begin() + old_sz);
        pool.count.push_back(1);

        // Compute light_side = min(popcount, n_tips - popcount)
        int pc = 0;
        for (int bin = 0; bin < n_bins; ++bin) {
          pc += count_bits(canon_buf[bin]);
        }
        pool.light_side.push_back(std::min(pc, n_tips - pc));

        // Insert pointer into map (points into pool.data)
        split_map[pool.split(idx)] = idx;
      }

      // Only record unique splits per tree
      if (members.empty() || members.back() != idx) {
        // Check if already present (small set, linear scan ok)
        bool found = false;
        for (int m : members) {
          if (m == idx) { found = true; break; }
        }
        if (!found) members.push_back(idx);
      }
    }
  }

  return pool;
}


// ============================================================================
// transfer_dist_mat: pairwise transfer distance between all unique splits
// ============================================================================

// Returns a flat n_splits x n_splits matrix (row-major).
// DIST[i * n + j] = transfer distance between split i and split j.
static std::vector<int> transfer_dist_mat(const PooledSplits& pool) {
  const int M = pool.n_splits;
  const int nb = pool.n_bins;
  const int nt = pool.n_tips;

  std::vector<int> dist(M * M, 0);

  for (int i = 0; i < M; ++i) {
    const splitbit* a = pool.split(i);
    for (int j = i + 1; j < M; ++j) {
      const splitbit* b = pool.split(j);
      int hamming = 0;
      for (int bin = 0; bin < nb; ++bin) {
        hamming += count_bits(a[bin] ^ b[bin]);
      }
      int td = std::min(hamming, nt - hamming);
      dist[i * M + j] = td;
      dist[j * M + i] = td;
    }
  }
  return dist;
}


// ============================================================================
// compute_td: transfer dissimilarity cost for each split
// ============================================================================

static std::vector<double> compute_td(
    const std::vector<int>& dist,
    const PooledSplits& pool,
    bool scale
) {
  const int M = pool.n_splits;
  const int n_tree = static_cast<int>(pool.tree_members.size());

  std::vector<double> td(M, 0.0);

  for (int t = 0; t < n_tree; ++t) {
    const auto& members = pool.tree_members[t];
    const int n_mem = static_cast<int>(members.size());

    for (int b = 0; b < M; ++b) {
      // Min distance from split b to any split in tree t
      int min_d = pool.light_side[b] - 1; // sentinel distance
      for (int k = 0; k < n_mem; ++k) {
        int d = dist[b * M + members[k]];
        if (d < min_d) min_d = d;
      }
      const int p_minus_1 = pool.light_side[b] - 1;
      if (p_minus_1 <= 0) continue;
      if (scale) {
        double contrib = static_cast<double>(min_d) / p_minus_1;
        td[b] += (contrib < 1.0) ? contrib : 1.0;
      } else {
        td[b] += (min_d < p_minus_1) ? min_d : p_minus_1;
      }
    }
  }
  return td;
}


// ============================================================================
// compat_mat: pairwise compatibility between all unique splits
// ============================================================================

// Returns flat M x M bool matrix.
// compat[i * M + j] = true iff splits i and j are compatible.
static std::vector<uint8_t> compat_mat(const PooledSplits& pool) {
  const int M = pool.n_splits;
  const int nb = pool.n_bins;
  const int last_bin = nb - 1;
  const splitbit lm = pool.last_mask;

  std::vector<uint8_t> compat(M * M, 1);

  for (int i = 0; i < M; ++i) {
    const splitbit* a = pool.split(i);
    for (int j = i + 1; j < M; ++j) {
      const splitbit* b = pool.split(j);
      bool ab = false, anb = false, nab = false, nanb = false;
      for (int bin = 0; bin < nb; ++bin) {
        splitbit mask = (bin == last_bin) ? lm : ~splitbit(0);
        splitbit a_bin = a[bin] & mask;
        splitbit b_bin = b[bin] & mask;
        if (!ab)   ab   = (a_bin & b_bin) != 0;
        if (!anb)  anb  = (a_bin & ~b_bin & mask) != 0;
        if (!nab)  nab  = (~a_bin & b_bin & mask) != 0;
        if (!nanb) nanb = (~a_bin & ~b_bin & mask) != 0;
        if (ab && anb && nab && nanb) break;
      }
      bool comp = !ab || !anb || !nab || !nanb;
      compat[i * M + j] = comp ? 1 : 0;
      compat[j * M + i] = comp ? 1 : 0;
    }
  }
  return compat;
}


// ============================================================================
// Greedy loop helpers
// ============================================================================

// Sentinel distance for split b (distance to a leaf bipartition)
inline int sent_dist(int b, const PooledSplits& pool) {
  return pool.light_side[b] - 1;
}

// "Effective distance" from split b to its match (or sentinel if match == -1)
inline int eff_dist(int b, int match_idx, const std::vector<int>& dist,
                    int M, const PooledSplits& pool) {
  return (match_idx < 0) ? sent_dist(b, pool) : dist[b * M + match_idx];
}

// Find second-closest included split to split b (excluding matchIdx).
// Returns {index, distance} as a pair; index = -1 means sentinel.
static std::pair<int, int> find_second(
    int b, int matchIdx,
    const std::vector<uint8_t>& incl,
    const std::vector<int>& dist,
    int M, const PooledSplits& pool,
    bool scale
) {
  int p_minus_1 = sent_dist(b, pool);
  int best = -1;
  int best_d = (p_minus_1 > 0) ? p_minus_1 : 0;

  for (int c = 0; c < M; ++c) {
    if (!incl[c] || c == matchIdx) continue;
    int d = dist[b * M + c];
    if (d < best_d) {
      best_d = d;
      best = c;
    }
  }
  if (best >= 0) {
    if (scale) {
      if (p_minus_1 <= 0 ||
          static_cast<double>(best_d) / p_minus_1 >= 1.0)
        return {-1, p_minus_1};
    } else {
      if (best_d >= p_minus_1)
        return {-1, p_minus_1};
    }
  }
  return {best, (best >= 0) ? best_d : p_minus_1};
}


// ============================================================================
// GreedyState: persistent state for the greedy loop.
//
// Maintains match_dist[], match2_dist[], n_incompat[], weight[] to avoid
// redundant O(M) scans in benefit calculations and compat checks.
// ============================================================================

struct GreedyState {
  int M;
  int n_tip;
  int n_incl;
  bool scale;

  std::vector<int>&       match;
  std::vector<int>&       match2;
  std::vector<uint8_t>&   incl;
  const std::vector<int>& dist;
  const std::vector<double>& td;
  const PooledSplits&     pool;
  const std::vector<uint8_t>& compat;

  // Cached per-split arrays
  std::vector<int>    match_dist;   // dist to match[b], or sentinel
  std::vector<int>    match2_dist;  // dist to match2[b], or sentinel
  std::vector<int>    n_incompat;   // # included splits incompatible with b
  std::vector<double> weight;       // count[b] / (light_side[b] - 1) or count[b]

  GreedyState(
      std::vector<int>& match_, std::vector<int>& match2_,
      std::vector<uint8_t>& incl_,
      const std::vector<int>& dist_, const std::vector<double>& td_,
      const PooledSplits& pool_, const std::vector<uint8_t>& compat_,
      bool scale_, int M_, int n_tip_
  ) : M(M_), n_tip(n_tip_), n_incl(0), scale(scale_),
      match(match_), match2(match2_), incl(incl_),
      dist(dist_), td(td_), pool(pool_), compat(compat_),
      match_dist(M_), match2_dist(M_), n_incompat(M_, 0), weight(M_)
  {
    // Precompute weight = count[b] / (light_side[b] - 1) for scaled,
    //                     count[b] for unscaled
    for (int b = 0; b < M; ++b) {
      int p1 = sent_dist(b, pool);
      if (p1 <= 0) { weight[b] = 0.0; continue; }
      weight[b] = scale
        ? static_cast<double>(pool.count[b]) / p1
        : static_cast<double>(pool.count[b]);
    }

    // Initialize match_dist / match2_dist from current match/match2
    for (int b = 0; b < M; ++b) {
      match_dist[b]  = eff_dist(b, match[b],  dist, M, pool);
      match2_dist[b] = eff_dist(b, match2[b], dist, M, pool);
    }

    // Initialize n_incompat from current incl
    for (int i = 0; i < M; ++i) {
      if (!incl[i]) continue;
      n_incl++;
      for (int j = 0; j < M; ++j) {
        if (!compat[j * M + i]) n_incompat[j]++;
      }
    }
  }

  // O(1) compat check
  bool is_compatible(int idx) const {
    return n_incompat[idx] == 0 && n_incl < n_tip - 3;
  }

  // Add-benefit using precomputed match_dist[] and weight[]
  double add_benefit(int c) const {
    double ben = -td[c];
    for (int b = 0; b < M; ++b) {
      int diff = match_dist[b] - dist[b * M + c];
      if (diff <= 0) continue;
      ben += diff * weight[b];
    }
    return ben;
  }

  // Remove-benefit using precomputed match2_dist[] and weight[]
  double remove_benefit(int c) const {
    double ben = td[c];
    for (int b = 0; b < M; ++b) {
      if (match[b] != c) continue;
      int diff = dist[b * M + c] - match2_dist[b];
      if (diff >= 0) continue;
      ben += diff * weight[b];
    }
    return ben;
  }

  // Execute add: update incl, match, match2, match_dist, match2_dist,
  //              n_incompat, n_incl
  void do_add(int idx) {
    incl[idx] = 1;
    n_incl++;

    // Update n_incompat for all splits
    for (int j = 0; j < M; ++j) {
      if (!compat[j * M + idx]) n_incompat[j]++;
    }

    // Update match/match2 and their cached distances
    for (int b = 0; b < M; ++b) {
      int new_d = dist[b * M + idx];
      if (new_d < match_dist[b]) {
        match2[b]      = match[b];
        match2_dist[b] = match_dist[b];
        match[b]       = idx;
        match_dist[b]  = new_d;
      } else if (new_d < match2_dist[b]) {
        match2[b]      = idx;
        match2_dist[b] = new_d;
      }
    }
  }

  // Execute remove: update incl, match, match2, match_dist, match2_dist,
  //                 n_incompat, n_incl
  void do_remove(int idx) {
    incl[idx] = 0;
    n_incl--;

    // Update n_incompat
    for (int j = 0; j < M; ++j) {
      if (!compat[j * M + idx]) n_incompat[j]--;
    }

    // Update match/match2 for affected splits
    for (int b = 0; b < M; ++b) {
      if (match[b] == idx) {
        // Promote match2 → match, find new match2
        match[b]      = match2[b];
        match_dist[b] = match2_dist[b];
        auto [sec, sec_d] = find_second(b, match[b], incl, dist, M, pool, scale);
        match2[b]      = sec;
        match2_dist[b] = sec_d;
      } else if (match2[b] == idx) {
        auto [sec, sec_d] = find_second(b, match[b], incl, dist, M, pool, scale);
        match2[b]      = sec;
        match2_dist[b] = sec_d;
      }
    }
  }
};


// ============================================================================
// Greedy "best" strategy
// ============================================================================

static void greedy_best(
    std::vector<int>& match,
    std::vector<int>& match2,
    std::vector<uint8_t>& incl,
    const std::vector<int>& dist,
    const std::vector<double>& td,
    const PooledSplits& pool,
    const std::vector<uint8_t>& compat,
    const std::vector<int>& sort_ord,
    bool scale,
    int M, int n_tip
) {
  GreedyState st(match, match2, incl, dist, td, pool, compat, scale, M, n_tip);

  while (true) {
    double best_ben = 0.0;
    int best_idx = -1;
    bool best_is_add = false;

    for (int si = 0; si < M; ++si) {
      int idx = sort_ord[si];
      if (incl[idx]) {
        double ben = st.remove_benefit(idx);
        if (ben > best_ben) {
          best_ben = ben;
          best_idx = idx;
          best_is_add = false;
        }
      } else {
        if (!st.is_compatible(idx)) continue;
        double ben = st.add_benefit(idx);
        if (ben > best_ben) {
          best_ben = ben;
          best_idx = idx;
          best_is_add = true;
        }
      }
    }

    if (best_ben <= 0.0 || best_idx < 0) break;

    if (best_is_add) {
      st.do_add(best_idx);
    } else {
      st.do_remove(best_idx);
    }
  }
}


// ============================================================================
// Greedy "first" strategy
// ============================================================================

static void greedy_first(
    std::vector<int>& match,
    std::vector<int>& match2,
    std::vector<uint8_t>& incl,
    const std::vector<int>& dist,
    const std::vector<double>& td,
    const PooledSplits& pool,
    const std::vector<uint8_t>& compat,
    const std::vector<int>& sort_ord,
    bool scale,
    int M, int n_tip
) {
  GreedyState st(match, match2, incl, dist, td, pool, compat, scale, M, n_tip);

  bool improving = true;
  while (improving) {
    improving = false;
    for (int si = 0; si < M; ++si) {
      int idx = sort_ord[si];
      if (incl[idx]) {
        if (st.remove_benefit(idx) > 0.0) {
          st.do_remove(idx);
          improving = true;
          break;
        }
      } else {
        if (!st.is_compatible(idx)) continue;
        if (st.add_benefit(idx) > 0.0) {
          st.do_add(idx);
          improving = true;
          break;
        }
      }
    }
  }
}


// ============================================================================
// Init matches (for init = "majority")
// ============================================================================

static void init_matches(
    std::vector<int>& match,
    std::vector<int>& match2,
    const std::vector<uint8_t>& incl,
    const std::vector<int>& dist,
    int M, const PooledSplits& pool,
    bool scale
) {
  // Collect included indices
  std::vector<int> inc_idx;
  for (int i = 0; i < M; ++i) {
    if (incl[i]) inc_idx.push_back(i);
  }
  if (inc_idx.empty()) return;

  for (int b = 0; b < M; ++b) {
    int p_minus_1 = pool.light_side[b] - 1;
    if (p_minus_1 <= 0) continue;

    // Find closest and second-closest among included
    int best = -1, second = -1;
    int best_d = p_minus_1, second_d = p_minus_1; // sentinel threshold

    for (int c : inc_idx) {
      int d = dist[b * M + c];
      if (d < best_d) {
        second = best;
        second_d = best_d;
        best = c;
        best_d = d;
      } else if (d < second_d) {
        second = c;
        second_d = d;
      }
    }

    if (best >= 0) {
      if (scale) {
        if (static_cast<double>(best_d) / p_minus_1 >= 1.0) continue;
      } else {
        if (best_d >= p_minus_1) continue;
      }
      match[b] = best;

      if (second >= 0) {
        if (scale) {
          if (static_cast<double>(second_d) / p_minus_1 < 1.0)
            match2[b] = second;
        } else {
          if (second_d < p_minus_1)
            match2[b] = second;
        }
      }
    }
  }
}


// ============================================================================
// Main exported function
// ============================================================================

//' Transfer consensus (C++ implementation)
//'
//' @param splits_list List of raw matrices (one per tree), each from as.Splits().
//' @param n_tip Number of tips.
//' @param scale Logical: use scaled transfer distance?
//' @param greedy_best Logical: TRUE for "best", FALSE for "first".
//' @param init_majority Logical: TRUE to start from majority-rule splits.
//'
//' @return A LogicalVector of length n_splits indicating which pooled splits
//'   are included in the consensus, plus attributes "raw_splits" (a raw matrix
//'   of all unique splits) and "light_side" (integer vector).
//' @keywords internal
// [[Rcpp::export]]
List cpp_transfer_consensus(
    const List& splits_list,
    const int n_tip,
    const bool scale,
    const bool greedy_best_flag,
    const bool init_majority
) {
  const int n_tree = splits_list.size();

  // ---- Pool unique splits ----
  PooledSplits pool = pool_splits(splits_list, n_tip);
  const int M = pool.n_splits;

  if (M == 0) {
    return List::create(
      Rcpp::Named("included") = LogicalVector(0),
      Rcpp::Named("raw_splits") = RawMatrix(0, 0),
      Rcpp::Named("light_side") = IntegerVector(0)
    );
  }

  // ---- Pairwise transfer distance matrix ----
  std::vector<int> dist = transfer_dist_mat(pool);

  // ---- Transfer dissimilarity cost ----
  std::vector<double> td = compute_td(dist, pool, scale);

  // ---- Compatibility matrix ----
  std::vector<uint8_t> compat = compat_mat(pool);

  // ---- Sort order (by count, descending) ----
  std::vector<int> sort_ord(M);
  std::iota(sort_ord.begin(), sort_ord.end(), 0);
  std::sort(sort_ord.begin(), sort_ord.end(),
            [&](int a, int b) { return pool.count[a] > pool.count[b]; });

  // ---- Mutable state ----
  std::vector<int> match(M, -1);   // -1 = sentinel
  std::vector<int> match2(M, -1);
  std::vector<uint8_t> incl(M, 0);

  // ---- Init from majority if requested ----
  if (init_majority) {
    double half = n_tree / 2.0;
    for (int i = 0; i < M; ++i) {
      if (pool.count[i] > half) {
        incl[i] = 1;
      }
    }
    init_matches(match, match2, incl, dist, M, pool, scale);
  }

  // ---- Greedy loop ----
  if (greedy_best_flag) {
    greedy_best(match, match2, incl, dist, td, pool, compat, sort_ord,
                scale, M, n_tip);
  } else {
    greedy_first(match, match2, incl, dist, td, pool, compat, sort_ord,
                 scale, M, n_tip);
  }

  // ---- Build output ----
  // Return raw split data for included splits so R can convert to phylo
  LogicalVector incl_r(M);
  for (int i = 0; i < M; ++i) incl_r[i] = incl[i] != 0;

  // Build raw matrix of ALL unique splits (n_splits x n_bytes)
  // The canonical form may differ from the original raw form: we need to
  // return the canonical splitbit data as raw bytes for as.phylo().
  const int n_bytes = pool.n_bins * static_cast<int>(sizeof(splitbit));
  RawMatrix raw_splits(M, n_bytes);
  for (int i = 0; i < M; ++i) {
    const unsigned char* src =
      reinterpret_cast<const unsigned char*>(pool.split(i));
    for (int j = 0; j < n_bytes; ++j) {
      raw_splits(i, j) = Rbyte(src[j]);
    }
  }

  IntegerVector light_side(M);
  for (int i = 0; i < M; ++i) light_side[i] = pool.light_side[i];

  return List::create(
    Rcpp::Named("included") = incl_r,
    Rcpp::Named("raw_splits") = raw_splits,
    Rcpp::Named("light_side") = light_side
  );
}


// ============================================================================
// Diagnostic: per-stage timing (for profiling without VTune)
// ============================================================================

//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericVector cpp_tc_profile(
    const List& splits_list,
    const int n_tip,
    const bool scale,
    const bool greedy_best_flag,
    const bool init_majority,
    const int n_iter
) {
  using clk = std::chrono::high_resolution_clock;
  const int n_tree = splits_list.size();

  double t_pool = 0, t_dist = 0, t_td = 0, t_compat = 0, t_greedy = 0;

  for (int iter = 0; iter < n_iter; ++iter) {
    auto t0 = clk::now();
    PooledSplits pool = pool_splits(splits_list, n_tip);
    auto t1 = clk::now();
    const int M = pool.n_splits;

    std::vector<int> dist = transfer_dist_mat(pool);
    auto t2 = clk::now();

    std::vector<double> td = compute_td(dist, pool, scale);
    auto t3 = clk::now();

    std::vector<uint8_t> compat_v = compat_mat(pool);
    auto t4 = clk::now();

    std::vector<int> sort_ord(M);
    std::iota(sort_ord.begin(), sort_ord.end(), 0);
    std::sort(sort_ord.begin(), sort_ord.end(),
              [&](int a, int b) { return pool.count[a] > pool.count[b]; });

    std::vector<int> match(M, -1);
    std::vector<int> match2(M, -1);
    std::vector<uint8_t> incl(M, 0);

    if (init_majority) {
      double half = n_tree / 2.0;
      for (int i = 0; i < M; ++i) {
        if (pool.count[i] > half) incl[i] = 1;
      }
      init_matches(match, match2, incl, dist, M, pool, scale);
    }

    if (greedy_best_flag) {
      greedy_best(match, match2, incl, dist, td, pool, compat_v, sort_ord,
                  scale, M, n_tip);
    } else {
      greedy_first(match, match2, incl, dist, td, pool, compat_v, sort_ord,
                   scale, M, n_tip);
    }
    auto t5 = clk::now();

    auto us = [](auto a, auto b) {
      return std::chrono::duration<double, std::micro>(b - a).count();
    };
    t_pool   += us(t0, t1);
    t_dist   += us(t1, t2);
    t_td     += us(t2, t3);
    t_compat += us(t3, t4);
    t_greedy += us(t4, t5);
  }

  double inv = 1.0 / n_iter;
  Rcpp::NumericVector result = {
    t_pool * inv, t_dist * inv, t_td * inv, t_compat * inv, t_greedy * inv
  };
  result.attr("names") = Rcpp::CharacterVector(
    {"pool_splits", "transfer_dist_mat", "compute_td", "compat_mat", "greedy"}
  );
  return result;
}
