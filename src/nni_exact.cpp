#include <Rcpp/Lightest>
#include <TreeTools/assert.h>

#include <algorithm> // for sort, unique, lower_bound
#include <array>
#include <numeric> // for iota
#include <unordered_set>

#include "nni_exact.h"

using namespace Rcpp;

namespace {

constexpr int32_t KEY_WORDS = 4;
constexpr int32_t MAX_SPLITS = NNI_EXACT_MAX_TIPS - 3;
static_assert(MAX_SPLITS <= KEY_WORDS * 4, "Key too small for splits");

using SplitArray = std::array<uint16_t, MAX_SPLITS>;

// A tree is identified by its sorted splits, packed four to a word.
struct TreeKey {
  uint64_t w[KEY_WORDS];
  bool operator==(const TreeKey& other) const {
    return w[0] == other.w[0] && w[1] == other.w[1] &&
      w[2] == other.w[2] && w[3] == other.w[3];
  }
};

struct TreeKeyHash {
  std::size_t operator()(const TreeKey& key) const {
    uint64_t h = 0x9E3779B97F4A7C15ULL;
    for (int32_t i = 0; i < KEY_WORDS; ++i) {
      h ^= key.w[i] + 0x9E3779B97F4A7C15ULL + (h << 6) + (h >> 2);
    }
    return static_cast<std::size_t>(h ^ (h >> 31));
  }
};

inline TreeKey pack(const SplitArray& s, const int32_t n_split) {
  TreeKey key{{0, 0, 0, 0}};
  for (int32_t i = 0; i < n_split; ++i) {
    key.w[i / 4] |= uint64_t(s[i]) << (16 * (i % 4));
  }
  return key;
}

inline void unpack(const TreeKey& key, const int32_t n_split, SplitArray& s) {
  for (int32_t i = 0; i < n_split; ++i) {
    s[i] = uint16_t(key.w[i / 4] >> (16 * (i % 4)));
  }
}

inline int32_t popcount16(uint16_t x) {
  int32_t n = 0;
  for (; x; x &= x - 1) ++n;
  return n;
}

// Replace split i of the sorted array s with `replacement`, keeping it sorted
inline TreeKey replaced(const SplitArray& s, const int32_t n_split,
                        const int32_t i, const uint16_t replacement) {
  SplitArray out = s;
  int32_t j = i;
  while (j > 0 && out[j - 1] > replacement) {
    out[j] = out[j - 1];
    --j;
  }
  while (j < n_split - 1 && out[j + 1] < replacement) {
    out[j] = out[j + 1];
    ++j;
  }
  out[j] = replacement;
  return pack(out, n_split);
}

// Each split S, polarized away from leaf 0, is a cluster of the tree rooted
// on leaf 0.  S has children X and Y; its parent P has children S and Q, and
// everything outside P (including leaf 0) forms the fourth subtree.
// The two NNIs across S's edge replace S with X | Q or with Y | Q.
template <typename Visit>
inline void for_each_neighbour(const SplitArray& s, const int32_t n_split,
                               const uint16_t all_but_0, Visit visit) {
  for (int32_t i = 0; i < n_split; ++i) {
    const uint16_t split = s[i];
    uint16_t child = 0, parent = all_but_0;
    int32_t child_size = 0, parent_size = popcount16(all_but_0);
    for (int32_t j = 0; j < n_split; ++j) {
      if (j == i) continue;
      const uint16_t other = s[j];
      if ((other & ~split) == 0) {
        const int32_t size = popcount16(other);
        if (size > child_size) {
          child = other;
          child_size = size;
        }
      } else if ((split & ~other) == 0) {
        const int32_t size = popcount16(other);
        if (size < parent_size) {
          parent = other;
          parent_size = size;
        }
      }
    }
    if (!child) {
      child = split & uint16_t(-split); // Lowest leaf
    }
    const uint16_t other_child = split ^ child;
    const uint16_t sibling = parent ^ split;
    visit(replaced(s, n_split, i, uint16_t(child | sibling)));
    visit(replaced(s, n_split, i, uint16_t(other_child | sibling)));
  }
}

} // namespace

int32_t nni_exact_small(const std::vector<uint16_t>& a,
                        const std::vector<uint16_t>& b,
                        const int32_t n_tip, const std::size_t budget) {
  ASSERT(n_tip <= NNI_EXACT_MAX_TIPS);
  const int32_t n_split = n_tip - 3;
  if (n_split < 1) return 0;
  ASSERT(int32_t(a.size()) == n_split && int32_t(b.size()) == n_split);

  SplitArray sa{}, sb{};
  std::copy(a.begin(), a.end(), sa.begin());
  std::copy(b.begin(), b.end(), sb.begin());
  std::sort(sa.begin(), sa.begin() + n_split);
  std::sort(sb.begin(), sb.begin() + n_split);
  const TreeKey key_a = pack(sa, n_split);
  const TreeKey key_b = pack(sb, n_split);
  if (key_a == key_b) return 0;

  const uint16_t all_but_0 = uint16_t(((1U << n_tip) - 1) ^ 1U);

  // Bidirectional breadth-first search.  The two balls stay disjoint, so the
  // first tree found in both lies on a shortest path.
  std::unordered_set<TreeKey, TreeKeyHash> seen_a{key_a}, seen_b{key_b};
  std::vector<TreeKey> front_a{key_a}, front_b{key_b}, next;
  int32_t depth_a = 0, depth_b = 0;
  SplitArray s{};

  while (true) {
    const bool from_a = front_a.size() <= front_b.size();
    std::vector<TreeKey>& front = from_a ? front_a : front_b;
    std::unordered_set<TreeKey, TreeKeyHash>& seen = from_a ? seen_a : seen_b;
    const std::unordered_set<TreeKey, TreeKeyHash>& other =
      from_a ? seen_b : seen_a;

    next.clear();
    bool met = false;
    bool over_budget = false;
    for (const TreeKey& key : front) {
      unpack(key, n_split, s);
      for_each_neighbour(s, n_split, all_but_0, [&](const TreeKey& nbr) {
        if (met || over_budget) return;
        if (other.count(nbr)) {
          met = true;
        } else if (seen.insert(nbr).second) {
          next.push_back(nbr);
          if (seen_a.size() + seen_b.size() > budget) {
            over_budget = true;
          }
        }
      });
      if (met) return depth_a + depth_b + 1;
      if (over_budget) return -1;
    }
    ASSERT(!next.empty()); // NNI tree space is connected
    front.swap(next);
    if (from_a) {
      ++depth_a;
    } else {
      ++depth_b;
    }
  }
}

namespace {

// Postorder edge matrix -> leaves below each node, n_bin words per node.
std::vector<uint64_t> node_splits(const IntegerMatrix& edge,
                                  const int32_t n_tip, const int32_t n_node,
                                  const int32_t n_bin) {
  std::vector<uint64_t> splits(std::size_t(n_node) * n_bin, 0);
  for (int32_t i = 0; i < n_tip; ++i) {
    splits[std::size_t(i) * n_bin + i / 64] = uint64_t(1) << (i % 64);
  }
  const int32_t n_edge = edge.nrow();
  for (int32_t i = 0; i < n_edge; ++i) {
    const std::size_t parent = std::size_t(edge(i, 0) - 1) * n_bin;
    const std::size_t child = std::size_t(edge(i, 1) - 1) * n_bin;
    for (int32_t j = 0; j < n_bin; ++j) {
      splits[parent + j] |= splits[child + j];
    }
  }
  return splits;
}

inline bool has_tip(const uint64_t* split, const int32_t tip) {
  return (split[tip / 64] >> (tip % 64)) & 1;
}

inline int32_t lowest_tip(const uint64_t* split, const int32_t n_bin) {
  for (int32_t j = 0; j < n_bin; ++j) {
    if (split[j]) {
      uint64_t x = split[j];
      int32_t bit = 0;
      while (!(x & 1)) {
        x >>= 1;
        ++bit;
      }
      return j * 64 + bit;
    }
  }
  return -1; // # nocov
}

inline int32_t lowest_absent_tip(const uint64_t* split, const int32_t n_tip) {
  for (int32_t i = 0; i < n_tip; ++i) {
    if (!has_tip(split, i)) return i;
  }
  return -1; // # nocov
}

// Unrooted view of a postorder edge matrix: a binary root is suppressed,
// its two children becoming adjacent.
struct UnrootedTree {
  int32_t n_tip, n_node, n_bin, root;
  bool rooted;
  std::vector<int32_t> parent;      // 0-based; -1 for the root
  std::vector<uint64_t> splits;

  UnrootedTree(const IntegerMatrix& edge, const int32_t n_tip_)
    : n_tip(n_tip_), n_node(edge.nrow() + 1), n_bin((n_tip_ + 63) / 64) {
    const int32_t n_edge = edge.nrow();
    root = edge(n_edge - 1, 0) - 1;
    parent.assign(n_node, -1);
    int32_t root_degree = 0;
    for (int32_t i = 0; i < n_edge; ++i) {
      parent[edge(i, 1) - 1] = edge(i, 0) - 1;
      if (edge(i, 0) - 1 == root) ++root_degree;
    }
    rooted = root_degree == 2;
    splits = node_splits(edge, n_tip, n_node, n_bin);
  }

  const uint64_t* split(const int32_t node) const {
    return &splits[std::size_t(node) * n_bin];
  }

  // Does the edge above internal `node` describe an internal edge of the
  // unrooted tree that no other node describes?  The two edges below a
  // binary root form one edge, described by `first_root_child`, and that
  // edge is internal only if both root children are internal.
  bool represents_edge(const int32_t node, const int32_t first_root_child)
    const {
    if (node == root) return false;
    if (rooted && parent[node] == root) {
      if (node != first_root_child) return false;
      for (int32_t i = n_tip; i < n_node; ++i) {
        if (i != node && parent[i] == root) return true;
      }
      return false;
    }
    return true;
  }
};

// Polarize so that tip 0 is absent
std::vector<uint64_t> normalized(const uint64_t* split, const int32_t n_tip,
                                 const int32_t n_bin) {
  std::vector<uint64_t> out(split, split + n_bin);
  if (out[0] & 1) {
    for (int32_t j = 0; j < n_bin; ++j) out[j] = ~out[j];
    const int32_t spare = n_bin * 64 - n_tip;
    if (spare) out[n_bin - 1] &= ~uint64_t(0) >> spare;
  }
  return out;
}

struct UnionFind {
  std::vector<int32_t> up;
  explicit UnionFind(const int32_t n) : up(n) {
    std::iota(up.begin(), up.end(), 0);
  }
  int32_t find(int32_t x) {
    while (up[x] != x) {
      up[x] = up[up[x]];
      x = up[x];
    }
    return x;
  }
  void unite(const int32_t a, const int32_t b) {
    up[find(a)] = find(b);
  }
};

} // namespace

int32_t nni_region_upper(const IntegerMatrix& edge1,
                         const IntegerMatrix& edge2,
                         const int32_t n_tip, const int32_t max_tips) {
  const UnrootedTree t1(edge1, n_tip);
  const UnrootedTree t2(edge2, n_tip);
  const int32_t n_node = t1.n_node;
  const int32_t n_bin = t1.n_bin;

  // The edge between the children of a binary root is represented by the
  // first child listed.
  int32_t first_child_1 = -1;
  for (int32_t i = 0; i < n_node; ++i) {
    if (t1.parent[i] == t1.root) { first_child_1 = i; break; }
  }
  int32_t first_child_2 = -1;
  for (int32_t i = 0; i < n_node; ++i) {
    if (t2.parent[i] == t2.root) { first_child_2 = i; break; }
  }

  // Splits of tree 2, sorted for lookup
  std::vector<std::vector<uint64_t>> splits2;
  for (int32_t node = n_tip; node < n_node; ++node) {
    if (t2.represents_edge(node, first_child_2)) {
      splits2.push_back(normalized(t2.split(node), n_tip, n_bin));
    }
  }
  std::sort(splits2.begin(), splits2.end());

  // Internal edges of tree 1 that do not occur in tree 2
  std::vector<bool> unmatched(n_node, false);
  for (int32_t node = n_tip; node < n_node; ++node) {
    if (t1.represents_edge(node, first_child_1)) {
      unmatched[node] = !std::binary_search(
        splits2.begin(), splits2.end(),
        normalized(t1.split(node), n_tip, n_bin));
    }
  }
  // Tree-2 splits absent from tree 1.  Shared splits would project
  // trivially onto every region; dropping them saves time.
  std::vector<std::vector<uint64_t>> matched1;
  for (int32_t node = n_tip; node < n_node; ++node) {
    if (t1.represents_edge(node, first_child_1) && !unmatched[node]) {
      matched1.push_back(normalized(t1.split(node), n_tip, n_bin));
    }
  }
  std::sort(matched1.begin(), matched1.end());
  std::vector<std::vector<uint64_t>> unmatched2;
  for (const auto& sp : splits2) {
    if (!std::binary_search(matched1.begin(), matched1.end(), sp)) {
      unmatched2.push_back(sp);
    }
  }
  
  // Neighbours of each internal node of tree 1, in the unrooted tree.
  // `node_of_edge` gives the node whose split describes the connecting edge.
  const int32_t sibling_of_first = [&]() {
    if (!t1.rooted) return -1;
    for (int32_t i = 0; i < n_node; ++i) {
      if (t1.parent[i] == t1.root && i != first_child_1) return i;
    }
    return -1; // # nocov
  }();
  std::vector<std::vector<int32_t>> nbr(n_node), nbr_edge(n_node);
  auto link = [&](const int32_t a, const int32_t b, const int32_t edge_node) {
    if (a >= n_tip) { nbr[a].push_back(b); nbr_edge[a].push_back(edge_node); }
    if (b >= n_tip) { nbr[b].push_back(a); nbr_edge[b].push_back(edge_node); }
  };
  for (int32_t node = 0; node < n_node; ++node) {
    const int32_t par = t1.parent[node];
    if (par < 0) continue;
    if (t1.rooted && par == t1.root) {
      if (node == first_child_1) link(node, sibling_of_first, node);
    } else {
      link(node, par, node);
    }
  }

  UnionFind regions(n_node);
  for (int32_t node = n_tip; node < n_node; ++node) {
    if (unmatched[node]) {
      for (std::size_t i = 0; i < nbr[node].size(); ++i) {
        if (nbr_edge[node][i] == node) regions.unite(node, nbr[node][i]);
      }
    }
  }
  std::vector<std::vector<int32_t>> members(n_node);
  for (int32_t node = n_tip; node < n_node; ++node) {
    if (t1.rooted && node == t1.root) continue;
    members[regions.find(node)].push_back(node);
  }

  int32_t total = 0;
  for (int32_t r = n_tip; r < n_node; ++r) {
    const std::vector<int32_t>& region = members[r];
    const int32_t n_edge = int32_t(region.size()) - 1;
    if (n_edge < 1) continue;
    const int32_t region_tips = n_edge + 3;
    if (n_edge == 1 || region_tips > max_tips) {
      total += nni_region_bound(n_edge); // Exact when n_edge == 1
      continue;
    }

    // One leaf to represent the subtree beyond each boundary edge
    std::vector<int32_t> rep;
    std::vector<int32_t> region_edges;
    for (const int32_t u : region) {
      for (std::size_t i = 0; i < nbr[u].size(); ++i) {
        const int32_t v = nbr[u][i];
        const int32_t edge_node = nbr_edge[u][i];
        if (edge_node >= n_tip && unmatched[edge_node]) {
          if (edge_node == u) region_edges.push_back(u);
          continue;
        }
        if (v < n_tip) {
          rep.push_back(v);
        } else if (t1.parent[v] == u) {
          rep.push_back(lowest_tip(t1.split(v), n_bin));
        } else {
          rep.push_back(lowest_absent_tip(t1.split(u), n_tip));
        }
      }
    }
    if (int32_t(rep.size()) != region_tips) {
      // Polytomy: fall back to the diameter
      total += nni_region_bound(n_edge);
      continue;
    }
    std::sort(rep.begin(), rep.end());
    ASSERT(int32_t(region_edges.size()) == n_edge);

    const uint16_t all = uint16_t((1U << region_tips) - 1);
    auto project = [&](const uint64_t* split,
                       std::vector<uint16_t>& out) {
      uint16_t mask = 0;
      for (int32_t j = 0; j < region_tips; ++j) {
        if (has_tip(split, rep[j])) mask |= uint16_t(1U << j);
      }
      if (mask & 1) mask ^= all;
      const int32_t size = popcount16(mask);
      if (size > 1 && size < region_tips - 1) out.push_back(mask);
    };
    std::vector<uint16_t> c0, c1;
    for (const int32_t node : region_edges) project(t1.split(node), c0);
    for (const auto& sp : unmatched2) project(sp.data(), c1);
    std::sort(c1.begin(), c1.end());
    c1.erase(std::unique(c1.begin(), c1.end()), c1.end());

    int32_t exact = -1;
    if (int32_t(c0.size()) == n_edge && int32_t(c1.size()) == n_edge) {
      exact = nni_exact_small(c0, c1, region_tips, NNI_EXACT_BUDGET);
    }
    total += exact < 0 ? nni_region_bound(n_edge) : exact;
  }
  return total;
}

// [[Rcpp::export]]
IntegerVector cpp_nni_exact(const IntegerMatrix& edge1,
                            const IntegerMatrix& edge2,
                            const IntegerVector& nTip,
                            const NumericVector& budget) {
  const int32_t n_tip = nTip[0];
  if (n_tip > NNI_EXACT_MAX_TIPS) {
    Rcpp::stop("Exact NNI distance is limited to %i leaves",
               NNI_EXACT_MAX_TIPS);
  }
  if (n_tip < 4) return IntegerVector::create(0);
  const UnrootedTree t1(edge1, n_tip);
  const UnrootedTree t2(edge2, n_tip);
  const uint16_t all = uint16_t((1U << n_tip) - 1);
  auto collect = [&](const UnrootedTree& t) {
    std::vector<uint16_t> out;
    for (int32_t node = n_tip; node < t.n_node; ++node) {
      if (node == t.root) continue;
      uint16_t mask = uint16_t(t.split(node)[0]);
      if (mask & 1) mask ^= all;
      const int32_t size = popcount16(mask);
      if (size > 1 && size < n_tip - 1) out.push_back(mask);
    }
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return out;
  };
  const std::vector<uint16_t> a = collect(t1), b = collect(t2);
  if (int32_t(a.size()) != n_tip - 3 || int32_t(b.size()) != n_tip - 3) {
    Rcpp::stop("Trees must be binary");
  }
  const int32_t dist = nni_exact_small(a, b, n_tip,
                                       std::size_t(budget[0]));
  return IntegerVector::create(dist < 0 ? NA_INTEGER : dist);
}
