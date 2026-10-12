#include <Rcpp/Lightest>
#include <algorithm> /* for std::sort */
#include <cmath> /* for std::log2 */
#include <vector>

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

// Minimum clustering-entropy contribution of a rooted subtree of each size;
// see .MinSubtreeEntropy(). Vectors are 1-indexed (element 0 unused).
static void min_subtree_entropy_impl(const int n, const std::vector<double> &h,
                                     std::vector<double> &cost,
                                     std::vector<int> &split) {
  cost.assign(n + 1, 0.0);
  split.assign(n + 1, 0);
  // rev[n - x] = cost[x], so that cost[i] + cost[s - i] reads both forwards
  std::vector<double> rev(n + 1, 0.0);
  if (n >= 2) {
    cost[2] = h[2];
    rev[n - 2] = cost[2];
    split[2] = 1;
  }
  for (int s = 3; s <= n; ++s) {
    const int half = s / 2;
    const double *a = &cost[1];          // a[k] = cost[k + 1]
    const double *b = &rev[n - s + 1];   // b[k] = cost[s - 1 - k]
    // Independent accumulators break the dependency chain; min is exact, so
    // the result does not depend on evaluation order
    double m0 = a[0] + b[0];
    double m1 = m0, m2 = m0, m3 = m0;
    int k = 1;
    for (; k + 3 < half; k += 4) {
      m0 = std::min(m0, a[k] + b[k]);
      m1 = std::min(m1, a[k + 1] + b[k + 1]);
      m2 = std::min(m2, a[k + 2] + b[k + 2]);
      m3 = std::min(m3, a[k + 3] + b[k + 3]);
    }
    for (; k < half; ++k) {
      m0 = std::min(m0, a[k] + b[k]);
    }
    const double bestDaughters = std::min(std::min(m0, m1), std::min(m2, m3));
    // Record the first (smallest) daughter attaining the minimum
    k = 0;
    while (a[k] + b[k] != bestDaughters) ++k;
    cost[s] = h[s] + bestDaughters;
    rev[n - s] = cost[s];
    split[s] = k + 1;
  }
}

// [[Rcpp::export]]
Rcpp::List min_subtree_entropy(const int n) {
  const std::vector<double> h = binary_entropies(n, xlog2x_table(n));
  std::vector<double> cost;
  std::vector<int> split;
  min_subtree_entropy_impl(n, h, cost, split);
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
  if (normalize) {
    min_subtree_entropy_impl(n, h, cost, split);
  }
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
    // cost[r] + cost[s] is bounded by its range minima, and by the cheapest
    // split of m leaves into two subtrees
    const double rsMin = std::max(costFloor[r0] + costFloor[s1],
                                  cost[m] - h[m]);
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
