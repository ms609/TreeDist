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
  if (n >= 2) {
    cost[2] = h[2];
    split[2] = 1;
  }
  for (int s = 3; s <= n; ++s) {
    const int half = s / 2;
    int best = 1;
    double bestDaughters = cost[1] + cost[s - 1];
    for (int i = 2; i <= half; ++i) {
      const double daughters = cost[i] + cost[s - i];
      if (daughters < bestDaughters) {
        bestDaughters = daughters;
        best = i;
      }
    }
    cost[s] = h[s] + bestDaughters;
    split[s] = best;
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

// Largest NNI clustering information distance for a single leaf count n >= 4;
// see .NNIMaxStep1().
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

  double best = -R_PosInf;
  int bestParts[4] = {0, 0, 0, 0};
  int bestSplits[2] = {0, 0};

  // Enumerate every partition of n into four parts p <= q <= r <= s.
  for (int p = 1; p <= n / 4; ++p) {
    for (int q = p; q <= (n - p) / 3; ++q) {
      const double h1 = h[p + q];
      const int rMax = (n - p - q) / 2;
      for (int r = q; r <= rMax; ++r) {
        const int s = n - p - q - r;
        const double h2 = h[p + r];
        const double h3 = h[p + s];
        // Drop the largest entropy: the NNI move keeps the two least balanced
        const double hMax = std::max(std::max(h1, h2), h3);
        const double keptEntropy = h1 + h2 + h3 - hMax;
        const int parts[4] = {p, q, r, s};
        const double joint = partition_entropy(parts, xlog2x, log2n, n);
        const double vi = 2 * joint - keptEntropy;
        const double value = normalize ?
          vi / (2 * (cost[p] + cost[q] + cost[r] + cost[s]) + keptEntropy) :
          vi;
        if (value > best) {
          best = value;
          bestParts[0] = p;
          bestParts[1] = q;
          bestParts[2] = r;
          bestParts[3] = s;
          // Discard the first split with the largest entropy
          const int sizes[3] = {p + q, p + r, p + s};
          const double entropies[3] = {h1, h2, h3};
          int drop = 0;
          for (int i = 1; i != 3; ++i) {
            if (entropies[i] > entropies[drop]) drop = i;
          }
          int k = 0;
          for (int i = 0; i != 3; ++i) {
            if (i != drop) bestSplits[k++] = sizes[i];
          }
          std::sort(bestSplits, bestSplits + 2);
        }
      }
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("value") = best,
    Rcpp::Named("subtrees") = Rcpp::IntegerVector(bestParts, bestParts + 4),
    Rcpp::Named("splits") = Rcpp::IntegerVector(bestSplits, bestSplits + 2)
  );
}
