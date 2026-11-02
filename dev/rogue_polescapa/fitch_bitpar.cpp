// Bit-parallel incremental Fitch: length added by each leaf of `tree` when it
// is coded `?` in every binary character. Prototype for a C++ cost estimate;
// see fitch_bench.R.
//
// Characters are packed 64 to a word as two bit planes (has-0, has-1).
// Downpass gives D[v] and the summed cost dc[v]; an uppass gives U[v], the
// Fitch set of the tree minus subtree(v) seen from v's attachment, and its
// summed cost uc[v]. Length with leaf x ambiguous = uc[x], so the added
// length is dc[root] - uc[x]. Cost: O(nNode * ceil(nChar / 64)).
#include <Rcpp.h>
#include <cstdint>
#include <vector>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector fitch_added(const IntegerMatrix edge, const int n_tip,
                          const LogicalMatrix chars) {
  // edge: preorder, rooted, binary; root = n_tip + 1. chars: nChar x nTip.
  const int n_edge = edge.nrow();
  const int n_node = n_edge + 1;
  const int n_char = chars.nrow();
  const int n_word = (n_char + 63) / 64;
  const int root = n_tip; // 0-based

  std::vector<int> parent(n_node, -1), kid1(n_node, -1), kid2(n_node, -1);
  for (int e = 0; e < n_edge; ++e) {
    const int p = edge(e, 0) - 1, c = edge(e, 1) - 1;
    parent[c] = p;
    (kid1[p] < 0 ? kid1[p] : kid2[p]) = c;
  }
  std::vector<uint64_t> d0(n_node * n_word, 0), d1(n_node * n_word, 0);
  std::vector<uint64_t> u0(n_node * n_word, 0), u1(n_node * n_word, 0);
  std::vector<int> dc(n_node, 0), uc(n_node, 0);
  std::vector<uint64_t> last(n_word, ~uint64_t(0));
  if (n_char % 64) last[n_word - 1] = (uint64_t(1) << (n_char % 64)) - 1;

  for (int t = 0; t < n_tip; ++t) {
    for (int c = 0; c < n_char; ++c) {
      const uint64_t bit = uint64_t(1) << (c % 64);
      (chars(c, t) ? d1 : d0)[t * n_word + c / 64] |= bit;
    }
  }
  for (int t = 0; t < n_tip; ++t) {
    for (int w = 0; w < n_word; ++w) {
      d0[t * n_word + w] = ~d1[t * n_word + w] & last[w];
    }
  }
  // Downpass: reverse preorder visits children before parents
  for (int e = n_edge - 1; e >= 0; --e) {
    const int v = edge(e, 0) - 1;
    if (edge(e, 1) - 1 != kid1[v]) continue; // once per internal node
    const int a = kid1[v], b = kid2[v];
    int cost = dc[a] + dc[b];
    for (int w = 0; w < n_word; ++w) {
      const uint64_t a0 = d0[a * n_word + w], a1 = d1[a * n_word + w],
        b0 = d0[b * n_word + w], b1 = d1[b * n_word + w];
      const uint64_t i0 = a0 & b0, i1 = a1 & b1;
      const uint64_t empty = ~(i0 | i1) & last[w];
      d0[v * n_word + w] = i0 | (empty & (a0 | b0));
      d1[v * n_word + w] = i1 | (empty & (a1 | b1));
      cost += __builtin_popcountll(empty);
    }
    dc[v] = cost;
  }
  // Uppass
  for (int e = 0; e < n_edge; ++e) {
    const int v = edge(e, 1) - 1, p = parent[v];
    const int s = kid1[p] == v ? kid2[p] : kid1[p];
    if (p == root) {
      for (int w = 0; w < n_word; ++w) {
        u0[v * n_word + w] = d0[s * n_word + w];
        u1[v * n_word + w] = d1[s * n_word + w];
      }
      uc[v] = dc[s];
      continue;
    }
    int cost = dc[s] + uc[p];
    for (int w = 0; w < n_word; ++w) {
      const uint64_t a0 = d0[s * n_word + w], a1 = d1[s * n_word + w],
        b0 = u0[p * n_word + w], b1 = u1[p * n_word + w];
      const uint64_t i0 = a0 & b0, i1 = a1 & b1;
      const uint64_t empty = ~(i0 | i1) & last[w];
      u0[v * n_word + w] = i0 | (empty & (a0 | b0));
      u1[v * n_word + w] = i1 | (empty & (a1 | b1));
      cost += __builtin_popcountll(empty);
    }
    uc[v] = cost;
  }
  IntegerVector added(n_tip);
  for (int t = 0; t < n_tip; ++t) added[t] = dc[root] - uc[t];
  return added;
}
