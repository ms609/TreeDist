#include <Rcpp/Lightest>
#include <algorithm> /* for std::sort */
#include <cmath> /* for std::log2 */
#include <vector>
#if defined(__SSE2__)
#include <emmintrin.h>
#elif defined(__aarch64__)
#include <arm_neon.h>
#endif

// x log2(x) for x = 0..n
inline std::vector<double> xlog2x_table(const int n) {
  std::vector<double> t(n + 1, 0.0);
  for (int x = 1; x <= n; ++x) {
    t[x] = x * std::log2(double(x));
  }
  return t;
}

// Entropy, in bits, of a partition of n = sum(parts) leaves.
// Arithmetic mirrors entropy_int() (information.h) so results are identical.
template <int N>
inline double partition_entropy(const int (&parts)[N],
                                const std::vector<double> &xlog2x,
                                const double log2n, const int n) {
  double sum = 0;
  for (int x : parts) {
    if (x > 0) {
      sum += xlog2x[x];
    }
  }
  return log2n - sum / n;
}

// Binary entropy of each k | (n - k) split, k = 0..n.
inline std::vector<double> binary_entropies(const int n,
                                            const std::vector<double> &xlog2x) {
  std::vector<double> h(n + 1, 0.0);
  if (n == 0) return h;
  const double log2n = std::log2(double(n));
  for (int k = 0; k <= n; ++k) {
    const int parts[2] = {k, n - k};
    h[k] = partition_entropy(parts, xlog2x, log2n, n);
  }
  return h;
}

// min(a[i] + b[i]) over 0 <= i < len, len >= 1.  min is exact, so the
// result does not depend on the order of evaluation.
inline double pair_sum_min(const double *a, const double *b, const int len) {
  int i = 0;
  double ret = a[0] + b[0];
#if defined(__SSE2__)
  if (len >= 8) {
    __m128d m0 = _mm_add_pd(_mm_loadu_pd(a), _mm_loadu_pd(b));
    __m128d m1 = _mm_add_pd(_mm_loadu_pd(a + 2), _mm_loadu_pd(b + 2));
    __m128d m2 = _mm_add_pd(_mm_loadu_pd(a + 4), _mm_loadu_pd(b + 4));
    __m128d m3 = _mm_add_pd(_mm_loadu_pd(a + 6), _mm_loadu_pd(b + 6));
    for (i = 8; i + 8 <= len; i += 8) {
      m0 = _mm_min_pd(m0, _mm_add_pd(_mm_loadu_pd(a + i), _mm_loadu_pd(b + i)));
      m1 = _mm_min_pd(m1, _mm_add_pd(_mm_loadu_pd(a + i + 2),
                                     _mm_loadu_pd(b + i + 2)));
      m2 = _mm_min_pd(m2, _mm_add_pd(_mm_loadu_pd(a + i + 4),
                                     _mm_loadu_pd(b + i + 4)));
      m3 = _mm_min_pd(m3, _mm_add_pd(_mm_loadu_pd(a + i + 6),
                                     _mm_loadu_pd(b + i + 6)));
    }
    const __m128d m = _mm_min_pd(_mm_min_pd(m0, m1), _mm_min_pd(m2, m3));
    ret = std::min(_mm_cvtsd_f64(m), _mm_cvtsd_f64(_mm_unpackhi_pd(m, m)));
  }
#elif defined(__aarch64__)
  if (len >= 8) {
    float64x2_t m0 = vaddq_f64(vld1q_f64(a), vld1q_f64(b));
    float64x2_t m1 = vaddq_f64(vld1q_f64(a + 2), vld1q_f64(b + 2));
    float64x2_t m2 = vaddq_f64(vld1q_f64(a + 4), vld1q_f64(b + 4));
    float64x2_t m3 = vaddq_f64(vld1q_f64(a + 6), vld1q_f64(b + 6));
    for (i = 8; i + 8 <= len; i += 8) {
      m0 = vminq_f64(m0, vaddq_f64(vld1q_f64(a + i), vld1q_f64(b + i)));
      m1 = vminq_f64(m1, vaddq_f64(vld1q_f64(a + i + 2), vld1q_f64(b + i + 2)));
      m2 = vminq_f64(m2, vaddq_f64(vld1q_f64(a + i + 4), vld1q_f64(b + i + 4)));
      m3 = vminq_f64(m3, vaddq_f64(vld1q_f64(a + i + 6), vld1q_f64(b + i + 6)));
    }
    ret = vminvq_f64(vminq_f64(vminq_f64(m0, m1), vminq_f64(m2, m3)));
  }
#endif
  // Independent accumulators break the dependency chain
  double m0 = ret, m1 = ret, m2 = ret, m3 = ret;
  for (; i + 4 <= len; i += 4) {
    m0 = std::min(m0, a[i] + b[i]);
    m1 = std::min(m1, a[i + 1] + b[i + 1]);
    m2 = std::min(m2, a[i + 2] + b[i + 2]);
    m3 = std::min(m3, a[i + 3] + b[i + 3]);
  }
  for (; i < len; ++i) {
    m0 = std::min(m0, a[i] + b[i]);
  }
  return std::min(std::min(m0, m1), std::min(m2, m3));
}

// First i minimizing a[i] + b[i] over 0 <= i < len, len >= 1; sets `min`.
// Minima of chunks locate the first minimum without a second full pass.
inline int pair_sum_argmin(const double *a, const double *b, const int len,
                           double &min) {
  constexpr int chunk = 64;
  if (len <= 16) {
    int at = 0;
    min = a[0] + b[0];
    for (int i = 1; i < len; ++i) {
      const double sum = a[i] + b[i];
      if (sum < min) {
        min = sum;
        at = i;
      }
    }
    return at;
  }
  min = R_PosInf;
  int start = 0;
  for (int c = 0; c < len; c += chunk) {
    const double chunkMin = pair_sum_min(a + c, b + c, std::min(chunk, len - c));
    if (chunkMin < min) {
      min = chunkMin;
      start = c;
    }
  }
  int i = start;
  while (a[i] + b[i] != min) ++i;
  return i;
}

// Lower bounds on cost[y] that are linear within aligned blocks of 2^j sizes,
// j >= jMin: the chord across the block, lowered by its largest excess over
// cost[].  A block is recorded once its last entry, or the last entry overall,
// is known.
class CostBounds {
  // (base, slope) of each block, levels concatenated; the bound on cost[y]
  // is base + slope * (y - block start)
  std::vector<double> table_;
  std::vector<int> offset_;
  int levels_;
  
public:
  static constexpr int jMin = 4;
  
private:
  void Record(const int j, const int a, const int e,
              const std::vector<double> &cost) {
    const double slope = e > a ? (cost[e] - cost[a]) / (e - a) : 0;
    double excess = 0;
    for (int y = a + 1; y < e; ++y) {
      excess = std::max(excess, cost[a] + slope * (y - a) - cost[y]);
    }
    double *entry = &table_[2 * (offset_[j] + (a >> j))];
    entry[0] = cost[a] - excess;
    entry[1] = slope;
  }
  
public:
  explicit CostBounds(const int n) : levels_(1) {
    while ((1 << levels_) <= n) ++levels_;
    offset_.assign(levels_ + 1, 0);
    for (int j = jMin; j < levels_; ++j) {
      offset_[j + 1] = offset_[j] + (n >> j) + 1;
    }
    table_.assign(2 * offset_[levels_], 0.0);
  }
  
  int levels() const { return levels_; }
  
  // cost[x] is now known; last: no further entries will follow
  void Complete(const int x, const std::vector<double> &cost, const bool last) {
    for (int j = jMin; j < levels_; ++j) {
      const int a = (x >> j) << j;
      if (x == a + (1 << j) - 1 || last) Record(j, a, x, cost);
    }
  }
  
  // Lower bound on min(cost[k] + cost[m - k]) over k0 <= k <= k1, using
  // level-j blocks, which must be recorded.  The bound is linear between
  // block boundaries, so is smallest at the end of a run.
  double PairMin(const int j, const int m, const int k0, const int k1) const {
    const double *t = &table_[2 * offset_[j]];
    auto Bound = [t, j](const int y) {
      const int b = y >> j;
      return t[2 * b] + t[2 * b + 1] * (y - (b << j));
    };
    double ret = R_PosInf;
    for (int k = k0; k <= k1; ) {
      const int runEnd = std::min(std::min(k1, k | ((1 << j) - 1)),
                                  m - (((m - k) >> j) << j));
      ret = std::min(ret, std::min(Bound(k) + Bound(m - k),
                                   Bound(runEnd) + Bound(m - runEnd)));
      k = runEnd + 1;
    }
    return ret;
  }
};

// Smaller daughter of the complete binary tree with s >= 2 leaves
inline int complete_split(const int s) {
  int p = 1; // largest power of two below s
  while (p * 2 < s) p *= 2;
  return std::max(p / 2, s - p);
}

// Minimum clustering-entropy contribution of a rooted subtree of each size;
// see .MinSubtreeEntropy(). Vectors are 1-indexed (element 0 unused).
//
// cost[s] = h[s] + min(cost[k] + cost[s - k]) over 1 <= k <= s / 2.  The
// first k attaining the minimum is recorded, as in an exhaustive scan.
// Two results proven in dev/min_subtree_entropy_proof.md avoid most of the
// scan:
//  - For s <= ceil(n / 2), the complete binary tree is the unique optimum, so
//    k = complete_split(s).
//  - For s >= 2n / 3, any k whose complete tree has a larger daughter of more
//    than m = n - s leaves is beaten by its own smaller daughter, leaving
//    k <= m + (largest power of two <= m).
// Remaining candidates are found by branch and bound over blocks of k, using
// lower bounds from `bounds`.  Returns whether `bounds` was filled; small n
// need none.
static bool min_subtree_entropy_impl(const int n, const std::vector<double> &h,
                                     std::vector<double> &cost,
                                     std::vector<int> &split,
                                     CostBounds &bounds) {
  constexpr double margin = 1e-9;
  // Visit blocks of 2^j sizes for j = jStart, jStart - 3, ..., starting with
  // 16-32 blocks; blocks that survive at the last level >= jMin are scanned.
  constexpr int jStep = 3;
  constexpr int jLeaf = CostBounds::jMin - 1;
  int jStart = -4; // floor(log2(s / 2)) - 4
  // Below this, an exhaustive scan is faster
  constexpr int minBlockHalf = 512;
  const bool track = n / 2 >= minBlockHalf;
  cost.assign(n + 1, 0.0);
  split.assign(n + 1, 0);
  // rev[n - x] = cost[x], so that cost[k] + cost[s - k] reads both forwards
  std::vector<double> rev(n + 1, 0.0);
  if (track) {
    bounds.Complete(0, cost, false);
    bounds.Complete(1, cost, false);
  }
  if (n >= 2) {
    cost[2] = h[2];
    rev[n - 2] = cost[2];
    split[2] = 1;
    if (track) bounds.Complete(2, cost, false);
  }
  for (int s = 3; s <= n; ++s) {
    if ((2 << (jStart + 4)) <= s / 2) ++jStart;
    if (2 * s <= n + 1) {
      const int k = complete_split(s);
      cost[s] = h[s] + (cost[k] + cost[s - k]);
      rev[n - s] = cost[s];
      split[s] = k;
      if (track) bounds.Complete(s, cost, s == n);
      continue;
    }
    int half = s / 2;
    if (3 * s >= 2 * n) {
      const int m = n - s;
      int q = 1;
      while (q * 2 <= m) q *= 2;
      half = std::min(half, m ? m + q : 1);
    }
    // The previous size's split is a good first guess
    int arg = std::min(std::max(split[s - 1], 1), half);
    double best = cost[arg] + cost[s - arg];
    if (arg < half && cost[arg + 1] + cost[s - arg - 1] < best) {
      ++arg;
      best = cost[arg] + cost[s - arg];
    }

    // Exhaustive scan of [k0, k1]
    auto Scan = [&](const int k0, const int k1) {
      const double *a = &cost[k0];        // a[i] = cost[k0 + i]
      const double *b = &rev[n - s + k0]; // b[i] = cost[s - k0 - i]
      double scanMin;
      const int i = pair_sum_argmin(a, b, k1 - k0 + 1, scanMin);
      if (scanMin < best || (scanMin == best && k0 + i < arg)) {
        best = scanMin;
        arg = k0 + i;
      }
    };
    auto Visit = [&](auto &self, const int j, const int k0, const int k1) {
      if (j <= jLeaf) {
        Scan(k0, k1);
        return;
      }
      // Bounds are valid for blocks that end before s
      const int known = ((s >> j) << j) - 1;
      for (int b0 = k0; b0 <= k1; b0 = (b0 | ((1 << j) - 1)) + 1) {
        const int b1 = std::min(b0 | ((1 << j) - 1), k1);
        if (s - b0 <= known && (b0 | ((1 << j) - 1)) <= known &&
            bounds.PairMin(j, s, b0, b1) > best + margin) {
          continue;
        }
        self(self, j - jStep, b0, b1);
      }
    };
    if (half < minBlockHalf) {
      Scan(1, half);
    } else {
      Visit(Visit, jStart, 1, half);
    }
    cost[s] = h[s] + best;
    rev[n - s] = cost[s];
    split[s] = arg;
    if (track) bounds.Complete(s, cost, s == n);
  }
  return track;
}

// [[Rcpp::export]]
Rcpp::List min_subtree_entropy(const int n) {
  const std::vector<double> h = binary_entropies(n, xlog2x_table(n));
  std::vector<double> cost;
  std::vector<int> split;
  CostBounds bounds(n);
  min_subtree_entropy_impl(n, h, cost, split, bounds);
  return Rcpp::List::create(
    Rcpp::Named("cost") = Rcpp::NumericVector(cost.begin() + 1, cost.end()),
    Rcpp::Named("split") = Rcpp::IntegerVector(split.begin() + 1, split.end())
  );
}

// Sum of the two smallest of three values
inline double two_smallest(const double a, const double b, const double c) {
  return a + b + c - std::max(std::max(a, b), c);
}

// Largest NNI clustering information distance for a single leaf count n >= 4;
// see .NNIMaxStep1().
//
// Candidates (p, q, r) are visited in the same lexicographic order as an
// exhaustive search, but a range of q, or an interval of r, is skipped
// wherever an upper bound on its values falls below the best value found so
// far (or a known lower bound on the optimum).  Skipped candidates cannot be
// optimal, so the result, including tie-breaking, is that of the exhaustive
// search.
//
// Bounds use an exact restatement of the variation of information.  Pairing
// p with y (and the remaining parts with each other) defines one of the three
// splits; its "deficit" D_y = (g(p, y) + g(rest)) / n, where
// g(a, b) = (a + b)(1 - Hb(a / (a + b))) is zero for a balanced pair and grows
// with imbalance.  Then VI = 2 - (sum of the two smallest deficits).
// g(a, b), a <= b, decreases in a and increases in b.
// [[Rcpp::export]]
Rcpp::List nni_max_step(const int n, const bool normalize) {
  const std::vector<double> xlog2x = xlog2x_table(n);
  const double log2n = std::log2(double(n));
  const std::vector<double> h = binary_entropies(n, xlog2x);
  std::vector<double> cost;
  std::vector<int> split;
  CostBounds bounds(normalize ? n : 0);
  const bool blockBounds = normalize &&
    min_subtree_entropy_impl(n, h, cost, split, bounds);
  // costFloor[x] = min(cost[x..n]): a lower bound on cost[y] for all y >= x.
  // pairFloor[m] likewise bounds the cheapest split of >= m leaves into two
  // subtrees, whose cost is cost[m] - h[m].
  std::vector<double> costFloor(cost);
  std::vector<double> pairFloor(cost.size());
  for (int x = int(cost.size()) - 1; x >= 0; --x) {
    pairFloor[x] = cost[x] - h[x];
    if (x + 1 < int(cost.size())) {
      costFloor[x] = std::min(costFloor[x], costFloor[x + 1]);
      pairFloor[x] = std::min(pairFloor[x], pairFloor[x + 1]);
    }
  }
  
  // n * deficit of the pair (a, b); 0 if a > b, i.e. outside its monotone range
  auto G = [&](const int a, const int b) {
    return a > b ? 0.0 : (a + b) - xlog2x[a + b] + xlog2x[a] + xlog2x[b];
  };
  
  auto Value = [&](const int p, const int q, const int r) {
    const int s = n - p - q - r;
    const double h1 = h[p + q];
    const double h2 = h[p + r];
    const double h3 = h[p + s];
    // Drop the largest entropy: the NNI move keeps the two least balanced
    const double hMax = std::max(std::max(h1, h2), h3);
    const double keptEntropy = h1 + h2 + h3 - hMax;
    const int parts[4] = {p, q, r, s};
    const double joint = partition_entropy(parts, xlog2x, log2n, n);
    const double vi = 2 * joint - keptEntropy;
    return normalize ?
      vi / (2 * (cost[p] + cost[q] + cost[r] + cost[s]) + keptEntropy) :
      vi;
  };
  
  // Upper bound on Value(p, q, r) for every q0 <= q <= q1, q <= r <= s.
  auto QBound = [&](const int p, const int q0, const int q1) {
    const int sMin = (n - p - q1 + 1) / 2;
    // r >= q0 and s >= sMin; g(q, s) is smallest at q = q1, s = sMin
    const double dq = G(p, q0);
    const double dr = std::max(G(p, q0), G(q1, sMin));
    const double ds = G(p, sMin);
    const double viMax = 2 - two_smallest(dq, dr, ds) / n;
    if (!normalize) return viMax;
    if (viMax <= 0) return 0.0;
    // Split sizes lie in [p + q0, n - 2 q0]; binary entropy is concave
    const double keptMin = 2 * std::min(h[p + q0], h[n - q0 - q0]);
    const double denom = 2 * (cost[p] + costFloor[q0] + pairFloor[n - p - q1]) +
      keptMin;
    return denom > 0 ? viMax / denom : R_PosInf;
  };
  
  // Upper bound on Value(p, q, r) for r0 <= r <= r1 <= (n - p - q) / 2.
  auto RBound = [&](const int p, const int q, const int r0, const int r1) {
    const int m = n - p - q;
    const int s0 = m - r0;
    const int s1 = m - r1;
    // Joint entropy rises as r approaches s
    const int parts[4] = {p, q, r1, s1};
    const double joint = partition_entropy(parts, xlog2x, log2n, n);
    // Binary entropy is concave, so is smallest at an interval's ends
    const double hMin = std::min(std::min(h[p + r0], h[p + r1]),
                                 std::min(h[p + s1], h[p + s0]));
    // The two kept entropies include at most one of h[p + q]
    const double keptMin = hMin + std::min(h[p + q], hMin);
    const double dq = std::max(G(p, q), G(r1, s1));
    const double dr = std::max(G(p, r0), G(q, s1));
    const double ds = std::max(G(p, s1), G(q, r0));
    const double viMax = std::min(2 * joint - keptMin,
                                  2 - two_smallest(dq, dr, ds) / n);
    if (!normalize) return viMax;
    if (viMax <= 0) return 0.0;
    // cost[r] + cost[s] is bounded by its range minima, by the cheapest
    // split of m leaves into two subtrees, and by block-wise bounds
    double rsMin = std::max(costFloor[r0] + costFloor[s1], cost[m] - h[m]);
    if (blockBounds) {
      int j = CostBounds::jMin;
      while ((1 << j) < r1 - r0 + 1 && j + 1 < bounds.levels()) ++j;
      rsMin = std::max(rsMin, bounds.PairMin(j, m, r0, r1));
    }
    const double denom = 2 * (cost[p] + cost[q] + rsMin) + keptMin;
    return denom > 0 ? viMax / denom : R_PosInf;
  };
  
  // Optima lie near p = q, r = s; scanning that line gives a lower bound on
  // the optimum.  The margin absorbs rounding in the bounds.
  constexpr double margin = 1e-10;
  double floorValue = -R_PosInf;
  for (int p = 1; p <= n / 4; ++p) {
    floorValue = std::max(floorValue, Value(p, p, (n - p - p) / 2));
  }
  floorValue -= margin;
  
  double best = -R_PosInf;
  int bestParts[4] = {0, 0, 0, 0};
  auto Threshold = [&]() { return std::max(best - margin, floorValue); };
  
  // Intervals are visited depth first, lower half first, so that candidates
  // are considered in lexicographic order.
  constexpr int leafSize = 8;
  int qStack[64][2];
  int rStack[64][2];
  
  // Enumerate every partition of n into four parts p <= q <= r <= s.
  for (int p = 1; p <= n / 4; ++p) {
    int qDepth = 0;
    qStack[0][0] = p;
    qStack[0][1] = (n - p) / 3;
    while (qDepth >= 0) {
      const int q0 = qStack[qDepth][0];
      const int q1 = qStack[qDepth][1];
      --qDepth;
      if (QBound(p, q0, q1) < Threshold()) continue;
      if (q0 != q1) {
        const int mid = q0 + (q1 - q0) / 2;
        ++qDepth;
        qStack[qDepth][0] = mid + 1;
        qStack[qDepth][1] = q1;
        ++qDepth;
        qStack[qDepth][0] = q0;
        qStack[qDepth][1] = mid;
        continue;
      }
      const int q = q0;
      int rDepth = 0;
      rStack[0][0] = q;
      rStack[0][1] = (n - p - q) / 2;
      while (rDepth >= 0) {
        const int r0 = rStack[rDepth][0];
        const int r1 = rStack[rDepth][1];
        --rDepth;
        if (RBound(p, q, r0, r1) < Threshold()) continue;
        if (r1 - r0 < leafSize) {
          for (int r = r0; r <= r1; ++r) {
            const double value = Value(p, q, r);
            if (value > best) {
              best = value;
              bestParts[0] = p;
              bestParts[1] = q;
              bestParts[2] = r;
              bestParts[3] = n - p - q - r;
            }
          }
          continue;
        }
        const int mid = r0 + (r1 - r0) / 2;
        ++rDepth;
        rStack[rDepth][0] = mid + 1;
        rStack[rDepth][1] = r1;
        ++rDepth;
        rStack[rDepth][0] = r0;
        rStack[rDepth][1] = mid;
      }
    }
  }
  
  // Report the two kept splits, discarding the first with the largest entropy
  const int p = bestParts[0];
  const int sizes[3] = {p + bestParts[1], p + bestParts[2], p + bestParts[3]};
  int drop = 0;
  for (int i = 1; i != 3; ++i) {
    if (h[sizes[i]] > h[sizes[drop]]) drop = i;
  }
  int bestSplits[2];
  int k = 0;
  for (int i = 0; i != 3; ++i) {
    if (i != drop) bestSplits[k++] = sizes[i];
  }
  std::sort(bestSplits, bestSplits + 2);
  
  return Rcpp::List::create(
    Rcpp::Named("value") = best,
    Rcpp::Named("subtrees") = Rcpp::IntegerVector(bestParts, bestParts + 4),
    Rcpp::Named("splits") = Rcpp::IntegerVector(bestSplits, bestSplits + 2)
  );
}
