#include <Rcpp/Lightest>
#include <TreeTools/assert.h> // for ASSERT
#include <TreeTools/types.h> // for intx
#include <TreeTools/keep_tip.h> // for keep_tip
#include <TreeTools/root_tree.h> // for root_on_node
#include <TreeTools/SplitList.h> // for SplitList
#include <TreeDist/popcount.h>
#include "reduce_tree.h"
#include "col_prices.h"
#include "split_counts.h"
#include <cmath> // for log2
#include <algorithm> // for sort, lower_bound, equal
#include <cstdint>
#include <vector>

using namespace Rcpp;

template <typename Pop>
IntegerVector calc_mismatch_size_impl(const RawMatrix& x, const RawMatrix& y) {
  
  const TreeTools::SplitList a(x), b(y);
  
  const int32 n_split = int32(x.rows());
  const int32 n_tip = x.attr("nTip");
  const int32 half_tip = n_tip / 2;
  const int32 last_bin = a.n_bins - 1;
  const int32 unset_tips = (n_tip % SL_BIN_SIZE)
    ? SL_BIN_SIZE - n_tip % SL_BIN_SIZE : 0;
  
  constexpr splitbit all_ones = ~(splitbit(0U));
  const splitbit unset_mask = all_ones >> unset_tips;
  
  IntegerVector ret(n_split * n_split);
  int *ret_ptr = ret.end();
  for (int32 bi = b.n_splits; bi--; ) {
    for (int32 ai = a.n_splits; ai--; ) {
      
      --ret_ptr;
      
      *ret_ptr = Pop::count(
        (a.state[ai][last_bin] ^ b.state[bi][last_bin]) & unset_mask
      );
      
      for (int32 bin = last_bin; bin--; ) {
        *ret_ptr += Pop::count(a.state[ai][bin] ^ b.state[bi][bin]);
      }
      
      if (*ret_ptr > half_tip) {
        *ret_ptr = n_tip - *ret_ptr;
      }
    }
  }
  return ret;
}

IntegerVector calc_mismatch_size(const RawMatrix& x, const RawMatrix& y) {
  return TreeDist::dispatch_popcount([&](auto pop) {
    return calc_mismatch_size_impl<decltype(pop)>(x, y);
  });
}

// [[Rcpp::export]]
IntegerVector mismatch_size (const RawMatrix& x, const RawMatrix& y) {
  ASSERT(x.rows() == y.rows() && "`x` and `y` differ in number of splits.");
  ASSERT(x.hasAttribute("nTip") && "`x` lacks nTip attribute");
  ASSERT(y.hasAttribute("nTip") && "`y` lacks nTip attribute");
  ASSERT(static_cast<int>(x.attr("nTip")) == static_cast<int>(y.attr("nTip"))
         && "`x` and `y` differ in `nTip`");
  return calc_mismatch_size(x, y);
}

template <typename Pop>
IntegerVector calc_confusion_impl(const RawMatrix &x, const RawMatrix &y) {
  
  const TreeTools::SplitList a(x), b(y);
  
  const int32 n_split = static_cast<int32>(x.rows());
  const int32 n_tip = x.attr("nTip");
  const int32 n_bin = a.n_bins;
  const int32 confusion_size = 4;
  
  IntegerVector ret(n_split * n_split * confusion_size);
  int *ret_ptr = ret.end();
  for (int32 bi = n_split; bi--; ) {
    const int32 nb = b.in_split[bi];
    const int32 nB = n_tip - nb;
    
    for (int32 ai = n_split; ai--; ) {
      
      // x divides tips into a|A; y divides tips into b|B
      int32 a_and_b = 0;
      for (int32 bin = n_bin; bin--; ) {
        a_and_b += Pop::count(a.state[ai][bin] & b.state[bi][bin]);
      }
      
      const int32 na = a.in_split[ai];
      const int32 a_and_B = na - a_and_b;
      const int32 A_and_b = nb - a_and_b;
      const int32 A_and_B = nB - a_and_B;
      
      *(--ret_ptr) = A_and_B;
      *(--ret_ptr) = A_and_b;
      *(--ret_ptr) = a_and_B;
      *(--ret_ptr) = a_and_b;
    }
  }
  ret.attr("dim") = Dimension(confusion_size, n_split, n_split);
  return ret;
}

IntegerVector calc_confusion(const RawMatrix &x, const RawMatrix &y) {
  return TreeDist::dispatch_popcount([&](auto pop) {
    return calc_confusion_impl<decltype(pop)>(x, y);
  });
}

// [[Rcpp::export]]
IntegerVector confusion(const RawMatrix& x, const RawMatrix& y) {
  ASSERT(x.rows() == y.rows() && "Input splits must contain same number of splits.");
  ASSERT(x.hasAttribute("nTip") && "`x` lacks nTip attribute");
  ASSERT(y.hasAttribute("nTip") && "`y` lacks nTip attribute");
  ASSERT(static_cast<int>(x.attr("nTip")) == static_cast<int>(y.attr("nTip"))
         && "`x` and `y` differ in `nTip`");
  return calc_confusion(x, y);
}

IntegerMatrix reverse (const IntegerMatrix x) {
  ASSERT(double(x.nrow()) <= double(std::numeric_limits<intx>::max())
         && "This many edges are not (yet) supported.");
  const intx n_edge = intx(x.nrow());
  ASSERT(n_edge % 2 == 0); // Tree is binary
  IntegerMatrix ret(n_edge, 2);
  
  for (intx i = n_edge / 2; i--; ) {
    const intx j = n_edge - i - 1;
    ret[i] = x[j];
    ret[i + n_edge] = x[j + n_edge];
    ret[j] = x[i];
    ret[j + n_edge] = x[i + n_edge];
  }
  return ret;
}

// tree1 and tree2 are binary trees in postorder with identical tip.labels
// [[Rcpp::export]]
List keep_and_reroot(const List& tree1,
                     const List& tree2,
                     const LogicalVector& keep) {
  
  IntegerMatrix postorder1 = tree1["edge"];
  ASSERT(postorder1.nrow() % 2 == 0); // Tree is binary
  IntegerMatrix pre1 = reverse(postorder1);
  
  IntegerMatrix postorder2 = tree2["edge"];
  ASSERT(postorder2.nrow() % 2 == 0); // Tree is binary
  IntegerMatrix pre2 = reverse(postorder2);
  
  ASSERT((postorder1.nrow() / 2 + 1) == keep.length());
  
  bool any_kept = false;
  for (auto i : keep) {
    if (i) {
      any_kept = true;
      break;
    }
  }
  
  if (!any_kept) {
    List nullTree = List::create(Named("edge") = IntegerMatrix(0, 2),
                                 _["Nnode"] = 0,
                                 _["tip.label"] = CharacterVector(0));
    
    nullTree.attr("class") = "phylo";
    nullTree.attr("order") = "preorder";
    return List::create(nullTree, nullTree);
  }
  
  IntegerMatrix ret_edge1 = TreeTools::keep_tip(pre1, keep);
  IntegerMatrix ret_edge2 = TreeTools::keep_tip(pre2, keep);
  
  const int32 n_edge = ret_edge1.nrow();
  const int32 n_node = n_edge / 2;
  
  if (n_node == 0) {
    const CharacterVector all_labels = tree1["tip.label"];
    const CharacterVector kept_labels = all_labels[keep];
    ASSERT(kept_labels.length() == 1);
    IntegerVector e_val = IntegerVector::create(2, 1);
    e_val.attr("dim") = Dimension(1, 2);
    List oneTipTree = List::create(
      Named("edge") = as<IntegerMatrix>(e_val),
      _["tip.label"] = kept_labels,
      _["Nnode"] = 1
    );
    
    oneTipTree.attr("class") = "phylo";
    oneTipTree.attr("order") = "preorder";
    return List::create(oneTipTree, oneTipTree);
  }
  
  const int32 n_tip = n_node + 1;
  CharacterVector
    old_label = tree1["tip.label"],
    new_labels(n_tip)
  ;
  
  int32 next_tip = n_tip;
  for (int32 i = int32(old_label.size()); i--; ) {
    if (keep[i]) {
      --next_tip;
      new_labels[next_tip] = old_label[i];
    }
  }
  
  List ret1 = List::create(Named("edge") = ret_edge1,
                           _["Nnode"] = n_node,
                           _["tip.label"] = new_labels),
       ret2 = List::create(Named("edge") = ret_edge2,
                           _["Nnode"] = n_node,
                           _["tip.label"] = new_labels);
  ret1.attr("class") = "phylo";
  ret1.attr("order") = "preorder";
  ret2.attr("class") = "phylo";
  ret2.attr("order") = "preorder";
  return List::create(
    TreeTools::root_on_node(ret1, 1),
    TreeTools::root_on_node(ret2, 1)
  );
}

// [[Rcpp::export]]
List keep_and_reduce(const List& tree1,
                     const List& tree2,
                     const LogicalVector& keep) {
  List rerooted = keep_and_reroot(tree1, tree2, keep);
  if (rerooted.size() == 1) {
    return Rcpp::List::create(R_NilValue);
  }
  
  List rerooted1 = rerooted[0];
  List rerooted2 = rerooted[1];
  IntegerMatrix edge1 = reverse(rerooted1["edge"]);
  IntegerMatrix edge2 = reverse(rerooted2["edge"]);
  
  if (edge1.nrow() < 1) {
    List nullTree = List::create(Named("edge") = NumericMatrix(0, 2),
                                 _["Nnode"] = 0,
                                 _["tip.label"] = CharacterVector(0));
    
    nullTree.attr("class") = "phylo";
    nullTree.attr("order") = "preorder";
    return List::create(nullTree, nullTree);
  }
  
  CharacterVector tip_label = rerooted1["tip.label"];
  return reduce_trees(edge1, edge2, tip_label);
}

// Leaf scoring for SPRDist(method = "Rogue").
//
// Follows the R implementation (.ScoreWithout() with the clustering
// information proxy) step by step; scores match to within numerical error.
// R-level overhead dominated its run time.

extern double mutual_clustering_bounded_score(const RawMatrix &x,
                                              const RawMatrix &y,
                                              const int32 n_tip,
                                              const double min_score,
                                              const ColPrices *prices);
extern double mutual_clustering_bounded_counts(const SplitCounts &sc,
                                               const double min_score,
                                               const ColPrices *prices);
extern std::vector<int_fast64_t> mutual_clustering_col_prices(
    const RawMatrix &x, const RawMatrix &y, const int32 n_tips);

namespace {

// .Machine$double.eps ^ 0.5
constexpr double SQRT_EPS = 1.4901161193847656e-08;
constexpr int SPLIT_BIN = 8;

struct TreeSplits {
  RawMatrix bits;
  double entropy; // As ClusteringEntropy()
};

using price_t = int_fast64_t;

// Prices (optimal column duals) for matching the splits of a pair of trees
// whose candidate leaves are being scored.  Each candidate's reduced trees
// are these trees restricted to fewer leaves, so each candidate split is the
// restriction of a split here, whose price is near-optimal for it.
struct ParentPrices {
  int n_tip = 0;
  int n_word = 0; // 64-bit words per split
  int n_split = 0;
  std::vector<uint64_t> bits; // Splits of tree 2, n_split * n_word
  std::vector<SEXP> label; // CHARSXP of each leaf
  RawMatrix splits1, splits2; // As as.Splits(), for computing prices

  // Prices are computed on first use: a parent LAP is wasted if the row
  // bound rejects every candidate.  Empty if there is nothing to match.
  const std::vector<price_t> &Price() const {
    if (!ready) {
      price = mutual_clustering_col_prices(splits1, splits2, n_tip);
      ready = true;
    }
    return price;
  }

private:
  mutable std::vector<price_t> price;
  mutable bool ready = false;
};

// Position of each of `labels` among the parent's leaves; false unless
// `labels` is a subsequence of the parent's labels.
bool MapLeaves(const ParentPrices &pp, const CharacterVector &labels,
               std::vector<int> &map) {
  const int n_tip = labels.size();
  map.resize(n_tip);
  int j = 0;
  for (int t = 0; t != n_tip; ++t) {
    const SEXP lab = STRING_ELT(labels, t);
    while (j < pp.n_tip && pp.label[j] != lab) ++j;
    if (j == pp.n_tip) return false;
    map[t] = j++;
  }
  return n_tip > 0;
}

// As TreeTools::as.Splits(): splits in node order, omitting the root and
// its non-tip child.  Edges must be in postorder.
// If `leaf_index` is given, also writes each split to `indexed` as
// `n_word` 64-bit words in which leaf t sets bit leaf_index[t].
// If `row_node` is given, also writes the node (1-based) of each split.
TreeSplits GetSplits(const IntegerMatrix &edge, const int n_tip,
                     const std::vector<int> *leaf_index = nullptr,
                     const int n_word = 0,
                     std::vector<uint64_t> *indexed = nullptr,
                     std::vector<int> *row_node = nullptr) {
  const int n_edge = edge.nrow();
  const int n_vertex = n_edge + 1;
  const int n_bin = (n_tip - 1) / SPLIT_BIN + 1;
  if (n_edge == n_tip || n_tip < 4) {
    return {RawMatrix(0, n_bin), 0};
  }
  const int *parent = edge.begin();
  const int *child = parent + n_edge;
  
  std::vector<uint8_t> below(n_vertex * n_bin, 0);
  for (int t = 0; t != n_tip; ++t) {
    below[t * n_bin + t / SPLIT_BIN] = uint8_t(1) << (t % SPLIT_BIN);
  }
  for (int e = 0; e != n_edge; ++e) {
    uint8_t *p = below.data() + (parent[e] - 1) * n_bin;
    const uint8_t *c = below.data() + (child[e] - 1) * n_bin;
    for (int w = 0; w != n_bin; ++w) {
      p[w] |= c[w];
    }
  }
  
  static thread_local std::vector<uint64_t> below_ix;
  if (leaf_index) {
    below_ix.assign(size_t(n_vertex) * n_word, 0);
    for (int t = 0; t != n_tip; ++t) {
      const int ix = (*leaf_index)[t];
      below_ix[size_t(t) * n_word + (ix >> 6)] = uint64_t(1) << (ix & 63);
    }
    for (int e = 0; e != n_edge; ++e) {
      uint64_t *p = below_ix.data() + size_t(parent[e] - 1) * n_word;
      const uint64_t *c = below_ix.data() + size_t(child[e] - 1) * n_word;
      for (int w = 0; w != n_word; ++w) {
        p[w] |= c[w];
      }
    }
  }
  
  const int root_node = parent[n_edge - 1];
  int root_child = child[n_edge - 1];
  int root_children = 1;
  for (int e = 0; e != n_edge - 1; ++e) {
    if (parent[e] == root_node) {
      ++root_children;
      if (child[e] > n_tip) {
        root_child = child[e];
      }
    }
  }
  const int trivial_origin = root_node - 1;
  const int trivial_two = root_children == 2 ? root_child - 1 : -1;
  
  std::vector<int> rows;
  rows.reserve(n_vertex - n_tip);
  for (int i = n_tip; i < n_vertex; ++i) {
    if (i != trivial_origin && i != trivial_two) {
      rows.push_back(i);
    }
  }
  const int n_split = rows.size();
  if (row_node) {
    row_node->resize(n_split);
    for (int r = 0; r != n_split; ++r) {
      (*row_node)[r] = rows[r] + 1;
    }
  }
  RawMatrix ret(n_split, n_bin);
  if (leaf_index) {
    indexed->resize(size_t(n_split) * n_word);
    for (int r = 0; r != n_split; ++r) {
      std::copy_n(below_ix.data() + size_t(rows[r]) * n_word, n_word,
                  indexed->data() + size_t(r) * n_word);
    }
  }
  // Sum of n * H(split), in bits
  const double n_lg2_n = n_tip * std::log2(double(n_tip));
  double n_entropy = 0;
  for (int r = 0; r != n_split; ++r) {
    const uint8_t *row = below.data() + rows[r] * n_bin;
    int count = 0;
    for (int j = 0; j != n_bin; ++j) {
      ret(r, j) = row[j];
      count += TreeDist::count_bits(row[j]);
    }
    const int other = n_tip - count;
    n_entropy += n_lg2_n - count * std::log2(double(count)) -
      other * std::log2(double(other));
  }
  return {ret, n_entropy / n_tip};
}

// Prices for the splits of tree 2 of a candidate, from those of its parent.
struct CandidatePrices {
  const ParentPrices *pp;
  const std::vector<int> *map; // Parent position of each candidate leaf
  const std::vector<uint64_t> *bits; // Candidate splits, parent leaf indexing
  int n_split;
  std::vector<price_t> price;
  
  // Each candidate split takes the least price among the parent splits that
  // restrict to it (a split of the candidate tree always has one, as the
  // candidate tree is the parent tree restricted to fewer leaves).
  static const price_t *Get(void *context) {
    CandidatePrices &cp = *static_cast<CandidatePrices *>(context);
    const ParentPrices &pp = *cp.pp;
    const std::vector<price_t> &parent_price = pp.Price();
    if (parent_price.empty()) return nullptr; // #nocov
    const std::vector<int> &map = *cp.map;
    const int n_word = pp.n_word;
    static thread_local std::vector<uint64_t> mask, keys;
    static thread_local std::vector<price_t> key_price;
    static thread_local std::vector<int> slot;
    mask.assign(n_word, 0);
    for (const int ix : map) {
      mask[ix >> 6] |= uint64_t(1) << (ix & 63);
    }
    // Polarize so that the candidate's first leaf lies outside each split
    const int ref_word = map[0] >> 6;
    const uint64_t ref_bit = uint64_t(1) << (map[0] & 63);
    int capacity = 4;
    while (capacity < 2 * pp.n_split) capacity <<= 1;
    slot.assign(capacity, -1);
    keys.resize(size_t(pp.n_split) * n_word);
    key_price.resize(pp.n_split);
    const auto Hash = [n_word, capacity](const uint64_t *k) {
      uint64_t h = 0x9E3779B97F4A7C15ULL;
      for (int w = 0; w != n_word; ++w) {
        h ^= k[w];
        h *= 0xBF58476D1CE4E5B9ULL;
        h ^= h >> 31;
      }
      return int(h & uint64_t(capacity - 1));
    };
    int n_key = 0;
    for (int q = 0; q != pp.n_split; ++q) {
      const uint64_t *src = &pp.bits[size_t(q) * n_word];
      uint64_t *k = &keys[size_t(n_key) * n_word];
      const uint64_t flip = (src[ref_word] & ref_bit) ? ~uint64_t(0) : 0;
      for (int w = 0; w != n_word; ++w) k[w] = (src[w] ^ flip) & mask[w];
      for (int h = Hash(k); ; h = (h + 1) & (capacity - 1)) {
        const int e = slot[h];
        if (e < 0) {
          slot[h] = n_key;
          key_price[n_key++] = parent_price[q];
          break;
        }
        if (std::equal(k, k + n_word, &keys[size_t(e) * n_word])) {
          key_price[e] = std::min(key_price[e], parent_price[q]);
          break;
        }
      }
    }
    cp.price.assign(cp.n_split, -1);
    std::vector<uint64_t> key(n_word);
    int n_found = 0;
    long double found_sum = 0;
    for (int i = 0; i != cp.n_split; ++i) {
      const uint64_t *src = cp.bits->data() + size_t(i) * n_word;
      const uint64_t flip = (src[ref_word] & ref_bit) ? ~uint64_t(0) : 0;
      for (int w = 0; w != n_word; ++w) key[w] = (src[w] ^ flip) & mask[w];
      for (int h = Hash(key.data()); slot[h] >= 0;
           h = (h + 1) & (capacity - 1)) {
        const int e = slot[h];
        if (std::equal(key.begin(), key.end(), &keys[size_t(e) * n_word])) {
          cp.price[i] = key_price[e];
          found_sum += key_price[e];
          ++n_found;
          break;
        }
      }
    }
    if (n_found == 0) return nullptr; // #nocov
    if (n_found < cp.n_split) {
      const price_t fill = price_t(found_sum / n_found); // #nocov
      for (price_t &v : cp.price) if (v < 0) v = fill; // #nocov
    }
    return cp.price.data();
  }
};

// As ClusteringInfoDist(tree1, tree2), for trees with identical tip labels
// and edges in postorder; but if the distance exceeds `cutoff` by more than
// numerical noise, may instead return a lower bound that exceeds `cutoff`.
// `pp` (may be null) holds prices for a parent pair of trees; `map` gives the
// parent position of each leaf, and `bits2` the splits of tree 2 in parent
// leaf indexing.
double BoundedCID(const TreeSplits &splits1, const TreeSplits &splits2,
                  const int n_tip, const double cutoff,
                  const ParentPrices *pp = nullptr,
                  const std::vector<int> *map = nullptr,
                  const std::vector<uint64_t> *bits2 = nullptr) {
  const double info = splits1.entropy + splits2.entropy;
  const double tol = std::max(1.0, info) * SQRT_EPS;
  double mci = 0;
  if (n_tip >= 4) {
    CandidatePrices cp{pp, map, bits2, int(splits2.bits.nrow()), {}};
    const ColPrices prices{&CandidatePrices::Get, &cp};
    mci = mutual_clustering_bounded_score(splits1.bits, splits2.bits, n_tip,
                                          (info - cutoff) / 2 - tol,
                                          pp ? &prices : nullptr);
  }
  const double ret = info - mci - mci;
  return ret < tol ? 0 : ret;
}

double BoundedCID(const List &tree1, const List &tree2, const double cutoff,
                  const ParentPrices *pp = nullptr) {
  const CharacterVector labels = tree1["tip.label"];
  const int n_tip = labels.size();
  std::vector<int> map;
  std::vector<uint64_t> bits2;
  if (pp && MapLeaves(*pp, labels, map)) {
    return BoundedCID(GetSplits(tree1["edge"], n_tip),
                      GetSplits(tree2["edge"], n_tip, &map, pp->n_word, &bits2),
                      n_tip, cutoff, pp, &map, &bits2);
  }
  return BoundedCID(GetSplits(tree1["edge"], n_tip),
                    GetSplits(tree2["edge"], n_tip), n_tip, cutoff);
}

// Index of the first split in `x` that also occurs in `table`, or -1;
// as TreeTools::FirstMatchingSplit().
int FirstMatchingSplit(const RawMatrix &x, const RawMatrix &table,
                       const int n_tip) {
  const int n_bin = x.ncol();
  const int n_spare = n_tip % SPLIT_BIN;
  const uint8_t last_mask = n_spare ? uint8_t((1u << n_spare) - 1u) : 0xFF;
  // Polarize each split so that leaf 1 lies outside it, then hash
  const auto Canonical = [&](const RawMatrix &splits, std::vector<uint8_t> &key,
                             std::vector<uint64_t> &hash) {
    const int n_row = splits.nrow();
    key.resize(n_row * n_bin);
    hash.resize(n_row);
    for (int r = 0; r != n_row; ++r) {
      const uint8_t flip = (splits(r, 0) & 1) ? 0xFF : 0;
      uint64_t h = 14695981039346656037ULL; // FNV-1a
      for (int j = 0; j != n_bin; ++j) {
        const uint8_t mask = j == n_bin - 1 ? last_mask : 0xFF;
        const uint8_t byte = (splits(r, j) ^ flip) & mask;
        key[r * n_bin + j] = byte;
        h = (h ^ byte) * 1099511628211ULL;
      }
      hash[r] = h;
    }
  };
  std::vector<uint8_t> x_key, t_key;
  std::vector<uint64_t> x_hash, t_hash;
  Canonical(x, x_key, x_hash);
  Canonical(table, t_key, t_hash);
  
  std::vector<std::pair<uint64_t, int>> sorted(t_hash.size());
  for (int k = 0; k != int(t_hash.size()); ++k) {
    sorted[k] = {t_hash[k], k};
  }
  std::sort(sorted.begin(), sorted.end());
  for (int i = 0; i != x.nrow(); ++i) {
    auto it = std::lower_bound(sorted.begin(), sorted.end(),
                               std::make_pair(x_hash[i], 0));
    for (; it != sorted.end() && it->first == x_hash[i]; ++it) {
      if (std::equal(x_key.begin() + i * n_bin,
                     x_key.begin() + (i + 1) * n_bin,
                     t_key.begin() + it->second * n_bin)) {
        return i;
      }
    }
  }
  return -1;
}

// As ReduceTrees(KeepTipPostorder(tree1, keep), KeepTipPostorder(tree2, keep));
// a list holding NULL if no tree remains after reduction.
List KeepAndReduce(const IntegerMatrix &pre1, const IntegerMatrix &pre2,
                   const CharacterVector &labels, const LogicalVector &keep) {
  int n_kept = 0;
  for (const int k : keep) {
    n_kept += k;
  }
  if (n_kept < 4) {
    // Unrooted trees on three leaves are identical, so reduce to nothing
    return List::create(R_NilValue);
  }
  return reduce_trees(
    reverse(TreeTools::root_binary(TreeTools::keep_tip(pre1, keep), 1)),
    reverse(TreeTools::root_binary(TreeTools::keep_tip(pre2, keep), 1)),
    labels[keep]);
}

// Trees with edges in postorder
List KeepAndReduce(const List &tree1, const List &tree2,
                   const LogicalVector &keep) {
  return KeepAndReduce(reverse(tree1["edge"]), reverse(tree2["edge"]),
                       tree1["tip.label"], keep);
}

// As .ScoreWithout(): the proxy distance between the trees that remain once
// the leaves not in `keep` are removed and the trees reduced.
// `pre1` and `pre2` are the trees' edges in preorder.
double ScoreWithout(const IntegerMatrix &pre1, const IntegerMatrix &pre2,
                    const CharacterVector &labels,
                    const LogicalVector &keep, const double cutoff,
                    const ParentPrices *pp) {
  const List outcome = KeepAndReduce(pre1, pre2, labels, keep);
  if (Rf_isNull(outcome[0])) {
    return R_NegInf;
  }
  const List out1 = outcome[0];
  const List out2 = outcome[1];
  const CharacterVector out_labels = out1["tip.label"];
  const int n_tip = out_labels.size();
  std::vector<int> map;
  std::vector<uint64_t> bits2;
  if (pp && !MapLeaves(*pp, out_labels, map)) {
    pp = nullptr; // #nocov
  }
  const TreeSplits splits1 = GetSplits(out1["edge"], n_tip);
  const TreeSplits splits2 = pp ?
    GetSplits(out2["edge"], n_tip, &map, pp->n_word, &bits2) :
    GetSplits(out2["edge"], n_tip);
  const int first_match = FirstMatchingSplit(splits1.bits, splits2.bits,
                                             n_tip);
  if (first_match < 0) {
    return BoundedCID(splits1, splits2, n_tip, cutoff, pp, &map, &bits2);
  }
  
  // Divide at the shared split, anchoring each half to the other
  LogicalVector half1(n_tip), half2(n_tip);
  for (int t = 0; t != n_tip; ++t) {
    half1[t] = (splits1.bits(first_match, t / SPLIT_BIN) >> (t % SPLIT_BIN)) & 1;
    half2[t] = !half1[t];
  }
  for (LogicalVector *half : {&half1, &half2}) {
    for (int t = 0; t != n_tip; ++t) {
      if (!(*half)[t]) {
        (*half)[t] = true;
        break;
      }
    }
  }
  const List sub1 = KeepAndReduce(out1, out2, half1);
  if (Rf_isNull(sub1[0])) {
    Rcpp::stop("Half-tree reduced to nothing; please report this bug."); // #nocov
  }
  const double score1 = BoundedCID(sub1[0], sub1[1], cutoff, pp);
  if (score1 > cutoff) {
    return score1;
  }
  const List sub2 = KeepAndReduce(out1, out2, half2);
  if (Rf_isNull(sub2[0])) {
    Rcpp::stop("Half-tree reduced to nothing; please report this bug."); // #nocov
  }
  return score1 + BoundedCID(sub2[0], sub2[1], cutoff - score1, pp);
}

// Scores candidates that drop a single leaf without rebuilding the trees.
//
// Removing leaf l from a binary tree suppresses its parent node, so exactly
// one split vanishes from each tree, and every other split restricts to a
// distinct split of the smaller tree.  If, once l is removed, the trees
// still cannot be reduced and still share no split, then the smaller trees
// are those that KeepAndReduce() would return (up to leaf numbering and
// rooting), so their mutual clustering cost matrix follows from the overlaps
// between the parent splits, without building, reducing or splitting trees.
class SingleDrop {
public:
  // `edge1` and `edge2`: postorder edges of the parent trees, binary and
  // rooted on leaf 1.  `bits1` and `bits2`: their splits, as from
  // GetSplits() with identity leaf indexing; `node1` and `node2`: the node
  // of each split.  Disabled if `n_tip` is zero.
  SingleDrop(const IntegerMatrix &edge1, const IntegerMatrix &edge2,
             const int n_tip, const int n_word,
             const std::vector<uint64_t> &bits1,
             const std::vector<uint64_t> &bits2,
             const std::vector<int> &node1, const std::vector<int> &node2)
    : n_tip_(n_tip), n_word_(n_word), n_split_(int(node1.size())),
      bits1_(bits1), bits2_(bits2) {
    ready_ = Init(edge1, edge2, node1, node2);
  }

  // Prepares to score the removal of `leaf` (0-based); false if it must be
  // scored by ScoreWithout() instead.
  bool Prepare(const int leaf) {
    const int tip = leaf + 1;
    if (!ready_ || tip == 1 || shared_[leaf] || Reduces(tip)) {
      return false;
    }
    const int drop1 = Vanishing(x_, row_of1_, tip, nullptr);
    drop2_ = Vanishing(y_, row_of2_, tip, &merge2_);
    const int n_child = n_split_ - 1;
    row_.resize(n_child);
    row_in_.resize(n_child);
    row_has_.resize(n_child);
    col_.resize(n_child);
    col_in_.resize(n_child);
    col_has_.resize(n_child);
    const int word = leaf >> 6;
    const uint64_t bit = uint64_t(1) << (leaf & 63);
    for (int r = 0, k = 0; r != n_split_; ++r) {
      if (r == drop1) continue;
      const int16_t has = (bits1_[size_t(r) * n_word_ + word] & bit) != 0;
      row_[k] = r;
      row_has_[k] = has;
      row_in_[k] = in1_[r] - has;
      ++k;
    }
    for (int c = 0, k = 0; c != n_split_; ++c) {
      if (c == drop2_) continue;
      const int16_t has = (bits2_[size_t(c) * n_word_ + word] & bit) != 0;
      col_[k] = c;
      col_has_[k] = has;
      col_in_[k] = in2_[c] - has;
      ++k;
    }
    return true;
  }

  // As ScoreWithout(), for the prepared candidate.
  double Score(const double cutoff, const ParentPrices *pp) {
    const int n_child = n_split_ - 1;
    const int n = n_tip_ - 1;
    // As GetSplits()
    const double n_lg2_n = n * std::log2(double(n));
    double entropy1 = 0;
    double entropy2 = 0;
    for (int k = 0; k != n_child; ++k) {
      const int a = row_in_[k];
      const int b = col_in_[k];
      entropy1 += n_lg2_n - a * std::log2(double(a)) -
        (n - a) * std::log2(double(n - a));
      entropy2 += n_lg2_n - b * std::log2(double(b)) -
        (n - b) * std::log2(double(n - b));
    }
    // As BoundedCID()
    const double info = entropy1 / n + entropy2 / n;
    const double tol = std::max(1.0, info) * SQRT_EPS;
    const SplitCounts sc{n, n_child, overlap_.data(), n_split_,
                         row_.data(), row_in_.data(), row_has_.data(),
                         col_.data(), col_in_.data(), col_has_.data()};
    pp_ = pp;
    const ColPrices prices{&SingleDrop::Prices, this};
    const double mci = mutual_clustering_bounded_counts(
      sc, (info - cutoff) / 2 - tol, pp ? &prices : nullptr);
    const double ret = info - mci - mci;
    return ret < tol ? 0 : ret;
  }

private:
  // False if no candidate can be scored here.
  bool Init(const IntegerMatrix &edge1, const IntegerMatrix &edge2,
            const std::vector<int> &node1, const std::vector<int> &node2) {
    const int n_tip = n_tip_;
    if (n_tip < 8 || n_split_ != n_tip - 3 || int(node2.size()) != n_split_) {
      return false;
    }
    const int n_vertex = 2 * n_tip;
    x_.Build(edge1, n_vertex);
    y_.Build(edge2, n_vertex);
    if (!x_.RootedOnLeaf1(n_tip) || !y_.RootedOnLeaf1(n_tip)) {
      return false;
    }
    row_of1_.assign(n_vertex, -1);
    row_of2_.assign(n_vertex, -1);
    in1_.resize(n_split_);
    in2_.resize(n_split_);
    for (int r = 0; r != n_split_; ++r) {
      row_of1_[node1[r]] = r;
      row_of2_[node2[r]] = r;
      in1_[r] = Count(&bits1_[size_t(r) * n_word_]);
      in2_[r] = Count(&bits2_[size_t(r) * n_word_]);
    }
    // Overlaps between parent splits; and, for each leaf, whether removing
    // it makes two splits identical, which happens if they differ in that
    // leaf alone
    overlap_.resize(size_t(n_split_) * n_split_);
    shared_.assign(n_tip, 0);
    const int spare = n_word_ * 64 - n_tip;
    const uint64_t last_mask = spare ? ~uint64_t(0) >> spare : ~uint64_t(0);
    for (int i = 0; i != n_split_; ++i) {
      const uint64_t *a = &bits1_[size_t(i) * n_word_];
      for (int j = 0; j != n_split_; ++j) {
        const uint64_t *b = &bits2_[size_t(j) * n_word_];
        int a_and_b = 0;
        for (int w = 0; w != n_word_; ++w) {
          a_and_b += TreeDist::count_bits(a[w] & b[w]);
        }
        overlap_[size_t(i) * n_split_ + j] = int16_t(a_and_b);
        const int a_xor_b = in1_[i] + in2_[j] - a_and_b - a_and_b;
        if (a_xor_b == 0 || a_xor_b == n_tip) {
          return false; // The trees share a split
        }
        if (a_xor_b == 1 || a_xor_b == n_tip - 1) {
          const uint64_t flip = a_xor_b == 1 ? 0 : ~uint64_t(0);
          for (int w = 0; w != n_word_; ++w) {
            const uint64_t diff = (a[w] ^ b[w] ^ flip) &
              (w == n_word_ - 1 ? last_mask : ~uint64_t(0));
            if (diff) {
              shared_[w * 64 + __builtin_ctzll(diff)] = 1;
              break;
            }
          }
        }
      }
    }
    return true;
  }

  // Tree structure, as in reduce_trees()
  struct Ledger {
    std::vector<int> child, sibling, parent;
    void Build(const IntegerMatrix &edge, const int n_vertex) {
      child.assign(n_vertex, 0);
      sibling.assign(n_vertex, 0);
      parent.assign(n_vertex, 0);
      for (int e = edge.nrow(); e--; ) {
        const int p = edge(e, 0);
        const int c = edge(e, 1);
        parent[c] = p;
        const int existing = child[p];
        if (existing) {
          sibling[existing] = c;
          sibling[c] = existing;
        } else {
          child[p] = c;
        }
      }
    }
    // Whether the root's children are leaf 1 and an internal node, as
    // reduce_trees() expects.
    bool RootedOnLeaf1(const int n_tip) const {
      const int root = parent[1];
      if (parent[root]) return false;
      const int other = child[root] == 1 ? sibling[1] : child[root];
      return sibling[1] == other && other > n_tip;
    }
    // As REMOVE_TIP() in reduce_trees(); `tip`'s parent is not the root.
    void RemoveTip(const int tip) {
      const int p = parent[tip];
      const int s = sibling[tip];
      parent[s] = parent[p];
      sibling[sibling[p]] = s;
      sibling[s] = sibling[p];
      child[parent[s]] = sibling[s];
      sibling[tip] = 0;
    }
  };

  int Count(const uint64_t *x) const {
    int ret = 0;
    for (int w = 0; w != n_word_; ++w) {
      ret += TreeDist::count_bits(x[w]);
    }
    return ret;
  }

  // The split that vanishes from `ledger`'s tree when `tip` is removed;
  // `merge`, if given, receives the split that it then duplicates, or -1.
  // `tip`'s parent is not the root.
  static int Vanishing(const Ledger &ledger, const std::vector<int> &row_of,
                       const int tip, int *merge) {
    const int parent = ledger.parent[tip];
    const int sibling = ledger.sibling[tip];
    if (row_of[parent] >= 0) {
      if (merge) *merge = row_of[sibling];
      return row_of[parent];
    }
    // `parent` is the root's child, whose split omits leaf 1 alone; the
    // sibling's split will omit only leaf 1 too
    if (merge) *merge = -1;
    return row_of[sibling];
  }

  // Whether reduce_trees() would remove further leaves once `tip` is gone:
  // the tests of its cherry, root and chain rules, applied to the parent
  // trees with `tip` removed.
  bool Reduces(const int tip) {
    xc_ = x_;
    yc_ = y_;
    xc_.RemoveTip(tip);
    yc_.RemoveTip(tip);
    const int n_tip = n_tip_;
    const auto IsTip = [n_tip](const int i) { return i <= n_tip; };
    const std::vector<int> &xs = xc_.sibling, &xp = xc_.parent,
      &xk = xc_.child, &ys = yc_.sibling, &yp = yc_.parent, &yk = yc_.child;
    // Cherries
    for (int i = 1; i <= n_tip; ++i) {
      const int s = xs[i];
      if (s && IsTip(s) && s == ys[i]) return true;
    }
    // Root
    int x_1 = xk[xs[1]];
    if (!IsTip(x_1)) x_1 = xs[x_1];
    if (IsTip(x_1)) {
      const int y_1 = yk[ys[1]];
      if (x_1 == y_1 || x_1 == ys[y_1]) return true;
    }
    // Chains
    for (int i = n_tip; i >= 2; --i) {
      if (!xs[i]) continue;
      const int aunt = xs[xp[i]];
      if (!aunt || !IsTip(aunt)) continue;
      const int g_aunt = xs[xp[aunt]];
      if (!g_aunt || !IsTip(g_aunt)) continue;
      const int gg_aunt = xs[xp[g_aunt]];
      if (aunt == ys[yp[i]] && g_aunt == ys[yp[aunt]]) {
        if (gg_aunt && IsTip(gg_aunt) && gg_aunt == ys[yp[g_aunt]]) {
          return true;
        }
      } else if (i == ys[yp[aunt]] && aunt == ys[yp[g_aunt]]) {
        const int niece1 = yk[ys[g_aunt]];
        if (gg_aunt && IsTip(gg_aunt) &&
            (gg_aunt == niece1 || gg_aunt == ys[niece1] ||
             gg_aunt == ys[g_aunt])) {
          return true;
        }
      }
    }
    return false;
  }

  // Parent prices for the prepared candidate's splits of tree 2.  A split
  // that two parent splits restrict to takes the lesser price.
  static const price_t *Prices(void *context) {
    SingleDrop &sd = *static_cast<SingleDrop *>(context);
    const std::vector<price_t> &parent = sd.pp_->Price();
    if (parent.empty()) return nullptr; // #nocov
    const int n_child = sd.n_split_ - 1;
    sd.price_.resize(n_child);
    for (int k = 0; k != n_child; ++k) {
      const int c = sd.col_[k];
      sd.price_[k] = c == sd.merge2_ ?
        std::min(parent[c], parent[sd.drop2_]) : parent[c];
    }
    return sd.price_.data();
  }

  const int n_tip_;
  const int n_word_;
  const int n_split_;
  const std::vector<uint64_t> &bits1_;
  const std::vector<uint64_t> &bits2_;
  bool ready_ = false;
  Ledger x_, y_, xc_, yc_;
  std::vector<int> row_of1_, row_of2_, in1_, in2_;
  std::vector<int16_t> overlap_;
  std::vector<uint8_t> shared_;
  // The prepared candidate
  std::vector<int> row_, row_in_, col_, col_in_;
  std::vector<int16_t> row_has_, col_has_;
  int drop2_ = -1;
  int merge2_ = -1;
  const ParentPrices *pp_ = nullptr;
  std::vector<price_t> price_;
};

}

// Scores for dropping each leaf from `tree1` and `tree2`, as .SPRRogue();
// `pair_scores` are scores for dropping each pair of leaves, in the order of
// combn(), if `depth` > 1 and every single-leaf score is finite.
// Candidates that cannot beat the best so far may receive a lower bound
// rather than their exact score.
// `incremental = FALSE` scores every candidate by rebuilding its trees
// (SingleDrop is not used), for testing.
// [[Rcpp::export]]
List spr_rogue_scores(const List &tree1, const List &tree2, const int depth,
                      const bool incremental = true) {
  const CharacterVector labels = tree1["tip.label"];
  const int n_tip = labels.size();
  const IntegerMatrix pre1 = reverse(tree1["edge"]);
  const IntegerMatrix pre2 = reverse(tree2["edge"]);
  
  // Prices from the full trees, to bound and warm-start each candidate
  ParentPrices pp;
  pp.n_tip = n_tip;
  pp.n_word = (n_tip + 63) / 64;
  pp.label.resize(n_tip);
  for (int t = 0; t != n_tip; ++t) {
    pp.label[t] = STRING_ELT(labels, t);
  }
  std::vector<int> identity(n_tip);
  for (int t = 0; t != n_tip; ++t) {
    identity[t] = t;
  }
  std::vector<uint64_t> bits1;
  std::vector<int> node1, node2;
  {
    const TreeSplits full1 = GetSplits(tree1["edge"], n_tip, &identity,
                                       pp.n_word, &bits1, &node1);
    const TreeSplits full2 = GetSplits(tree2["edge"], n_tip, &identity,
                                       pp.n_word, &pp.bits, &node2);
    pp.splits1 = full1.bits;
    pp.splits2 = full2.bits;
    pp.n_split = full2.bits.nrow();
  }
  const ParentPrices *pp_ptr = pp.n_split ? &pp : nullptr;
  SingleDrop single(tree1["edge"], tree2["edge"], incremental ? n_tip : 0,
                    pp.n_word, bits1, pp.bits, node1, node2);
  NumericVector scores(n_tip);
  double best = R_PosInf;
  bool all_finite = true;
  LogicalVector keep(n_tip, true);
  for (int i = 0; i != n_tip; ++i) {
    Rcpp::checkUserInterrupt();
    if (single.Prepare(i)) {
      scores[i] = single.Score(best, pp_ptr);
    } else {
      keep[i] = false;
      scores[i] = ScoreWithout(pre1, pre2, labels, keep, best, pp_ptr);
      keep[i] = true;
    }
    if (!R_FINITE(scores[i])) {
      all_finite = false;
      break;
    }
    best = std::min(best, double(scores[i]));
  }
  
  if (depth < 2 || !all_finite) {
    return List::create(Named("scores") = scores,
                        _["pair_scores"] = R_NilValue);
  }
  
  NumericVector pair_scores(n_tip * (n_tip - 1) / 2);
  int pair = 0;
  for (int i = 0; i != n_tip - 1; ++i) {
    Rcpp::checkUserInterrupt();
    keep[i] = false;
    for (int j = i + 1; j != n_tip; ++j) {
      keep[j] = false;
      pair_scores[pair] = ScoreWithout(pre1, pre2, labels, keep, best,
                                       pp_ptr);
      keep[j] = true;
      if (!R_FINITE(pair_scores[pair])) {
        return List::create(Named("scores") = scores,
                            _["pair_scores"] = pair_scores);
      }
      best = std::min(best, double(pair_scores[pair]));
      ++pair;
    }
    keep[i] = true;
  }
  return List::create(Named("scores") = scores,
                      _["pair_scores"] = pair_scores);
}

// Test hook for BoundedCID().
// [[Rcpp::export]]
double spr_bounded_cid(const List &tree1, const List &tree2,
                       const double cutoff) {
  return BoundedCID(tree1, tree2, cutoff);
}

// Leaves on one side of the first split in `tree1` that also occurs in
// `tree2`, as as.logical(as.Splits(tree1)[[FirstMatchingSplit(...)]]);
// a zero-length vector if the trees share no split.
// Trees must have identical tip labels and edges in postorder.
// [[Rcpp::export]]
LogicalVector spr_first_shared_split(const List &tree1, const List &tree2) {
  const CharacterVector labels = tree1["tip.label"];
  const int n_tip = labels.size();
  const RawMatrix splits1 = GetSplits(tree1["edge"], n_tip).bits;
  const int first_match = FirstMatchingSplit(
    splits1, GetSplits(tree2["edge"], n_tip).bits, n_tip);
  if (first_match < 0) {
    return LogicalVector(0);
  }
  LogicalVector ret(n_tip);
  for (int t = 0; t != n_tip; ++t) {
    ret[t] = (splits1(first_match, t / SPLIT_BIN) >> (t % SPLIT_BIN)) & 1;
  }
  return ret;
}

// As keep_and_reduce(), via a leaner route; returns a list whose first
// element is NULL if no tree remains after reduction.
// [[Rcpp::export]]
List spr_keep_and_reduce(const List &tree1, const List &tree2,
                         const LogicalVector &keep) {
  return KeepAndReduce(tree1, tree2, keep);
}
