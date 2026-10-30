#include <Rcpp/Lightest>
#include <TreeTools/SplitList.h>
#include <TreeTools/assert.h>

#include <algorithm> // for max
#include <cmath>
#include <memory> // for unique_ptr, uninitialized_fill_n
#include <mutex> // for once_flag
#include <vector>

#include "tree_distances.h"
#include "li_diameters.h"

using namespace Rcpp;
using TreeTools::SplitList;

constexpr int32_t NA_INT32 = -9999; //std::numeric_limits<int32_t>::min();


template <typename T, std::size_t StackSize>
class HybridBuffer {
public:
  explicit HybridBuffer(std::size_t n) : n_(n), data_(nullptr) {
    if (n_ <= StackSize) {
      data_ = stack_;
      std::uninitialized_fill_n(data_, n_, T{});
    } else {
      heap_ = std::make_unique<T[]>(n_);                               // #nocov
      data_ = heap_.get();                                             // #nocov
    }
  }

  ~HybridBuffer() {
    if (n_ <= StackSize) {
      std::destroy_n(data_, n_);
    }
  }

  // bounds-checked accessors (ASSERT is a macro in your project)
  T& operator[](std::ptrdiff_t i) {
    ASSERT(i >= 0 && static_cast<std::size_t>(i) < n_);
    return data_[i];
  }
  const T& operator[](std::ptrdiff_t i) const {
    ASSERT(i >= 0 && static_cast<std::size_t>(i) < n_);
    return data_[i];
  }

  T* data() {
    return data_;
  }
  const T* data() const {
    return data_;
  }
  std::size_t size() const {
    return n_;
  }

  bool on_stack() const { return n_ <= StackSize; }

private:
  std::size_t n_;
  T* data_;
  T stack_[StackSize];
  std::unique_ptr<T[]> heap_;
};


constexpr int32_t NNI_STACK_BINS = 16;
constexpr int32_t NNI_STACK_SPLITS = 512;
constexpr int32_t NNI_STACK_TIPS = NNI_STACK_BINS * NNI_STACK_SPLITS;

constexpr int32_t NNI_MAX_TIPS = 32768;
// If updating NNI_MAX_TIPS, also update lg2_ceiling constructor
// Once we get near int16_max, we need to worry about upgrading to int64
// where we're looking at products of ntips and nbins.

int32_t lg2_ceiling[NNI_MAX_TIPS + 1];
int32_t fack_lookup[NNI_MAX_TIPS + 1];
int32_t li[NNI_MAX_TIPS + 1];

static std::once_flag cache_init_flag;
static void initialize_cache() {
  lg2_ceiling[0] = -1;
  lg2_ceiling[1] = 0;
  lg2_ceiling[2] = 1;
  lg2_ceiling[3] = 2;
  lg2_ceiling[4] = 2;
  
  for (int32_t i = 4    + 1; i != 8    + 1; i++) lg2_ceiling[i] = 3;
  for (int32_t i = 8    + 1; i != 16   + 1; i++) lg2_ceiling[i] = 4;
  for (int32_t i = 16   + 1; i != 32   + 1; i++) lg2_ceiling[i] = 5;
  for (int32_t i = 32   + 1; i != 64   + 1; i++) lg2_ceiling[i] = 6;
  for (int32_t i = 64   + 1; i != 128  + 1; i++) lg2_ceiling[i] = 7;
  for (int32_t i = 128  + 1; i != 256  + 1; i++) lg2_ceiling[i] = 8;
  for (int32_t i = 256  + 1; i != 512  + 1; i++) lg2_ceiling[i] = 9;
  for (int32_t i = 512  + 1; i != 1024 + 1; i++) lg2_ceiling[i] = 10;
  for (int32_t i = 1024 + 1; i != 2048 + 1; i++) lg2_ceiling[i] = 11;
  for (int32_t i = 2048 + 1; i != 4096 + 1; i++) lg2_ceiling[i] = 12;
  for (int32_t i = 4096 + 1; i != 8192 + 1; i++) lg2_ceiling[i] = 13;
  for (int32_t i = 8192 + 1; i != 16384 + 1; i++) lg2_ceiling[i] = 14;
  for (int32_t i = 16384 + 1; i != 32768 + 1; i++) lg2_ceiling[i] = 15;
  
  /* Fack et al. (2002), theorem 3: with N = n - 2 internal nodes,
   * diam <= N ceil(lg N) - 2^ceil(lg N) + 1 + 2 (N - ceil(lg(N + 1))).
   * Never looser than the bound in their theorem 1. */
  for (int32_t i = 4; i != NNI_MAX_TIPS + 1; i++) {
    const int32_t n_node = i - 2;
    const int32_t lg_n = lg2_ceiling[n_node];
    fack_lookup[i] = n_node * lg_n - (int32_t(1) << lg_n) + 1 +
      2 * (n_node - lg2_ceiling[n_node + 1]);
  }
  for (int32_t i = 4; i != NNI_MAX_TIPS + 1; i++) {
    const int32_t log_ceiling = lg2_ceiling[i];
    const int32_t sorting_number = int32_t(
      i * log_ceiling - std::pow(2, log_ceiling) + 1);
  
    /* Calculate the shortest length of the longest path in a tree.
     * To make this path as short as possible, divide tips into three 
     * balanced trees, joined by a single node that will form part of every 
     * longest path.  One of these subtrees will be filled with >= n/3 nodes */
    const int32_t nodes_in_full = int32_t(std::ceil(log2(double(i) / 3)));
    /* We want to put a power of two tips in this subtree, such that every node is 
     * equally close to its root */
    const int32_t tips_in_full = int32_t(std::pow(2, nodes_in_full));
    /* Now the remaining tips must be spread sub-evenly between the remaining 
     * edges from this node.  Picture halving the tips; removing tips from one side
     * until it is a power of two will reduce the number of nodes by one, whilst 
     * at worst (if this brings the other side over a power of two) increasing the 
     * other side by one. */
    const int32_t tips_left = i - tips_in_full;
    // (log2(tips_left / 2) + 1) == (log2(tips_left) - 1) + 1
    const int32_t min_backbone_nodes = nodes_in_full + lg2_ceiling[tips_left];
    
    /* The worst-case scenario requires a move for every node not on the backbone: */
    const int32_t n_node = i - 2;
    const int32_t degenerate_distance = n_node - min_backbone_nodes;
  
    li[i] = sorting_number + degenerate_distance + degenerate_distance;
  }
}

inline void ensure_cache() {
  std::call_once(cache_init_flag, initialize_cache);
}

// Score subtree, add to score, and reset subtree size
void update_score(const int32_t subtree_edges, int32_t *lower_bound,
                  int32_t *tight_bound, 
                  int32_t *closest_upper_bound,
                  int32_t *loose_bound, 
                  int32_t *li_bound, int32_t *fack_bound) {
  if (subtree_edges) {
    const int32_t subtree_tips = subtree_edges + 3;
    // Each unmatched edge must undergo at least one NNI (Li et al. 1996).
    // Diameters bound a single pair's distance from above only, so they
    // cannot tighten this lower bound.
    const int32_t lower = subtree_edges;
    const int32_t this_li = li[subtree_tips],
      fack = fack_lookup[subtree_tips],
      upper = this_li < fack ? this_li : fack;
    
    *lower_bound += lower;
    *li_bound += this_li;
    *fack_bound += fack;
    *loose_bound += upper;
    
    if (subtree_tips <= N_EXACT) {
      const int32_t exact = exact_diameter[subtree_tips];
      if (*tight_bound != NA_INT32) {
        *tight_bound += exact;
      }
      *closest_upper_bound += exact;
    } else {
      *tight_bound = NA_INT32;
      *closest_upper_bound += upper;
    }
  }
}


constexpr int32_t NOT_TRIVIAL = std::numeric_limits<int32_t>::max();
constexpr int32_t N_NNI_BOUNDS = 6;

// Edges in postorder, with 1-based node numbers
struct NNITree {
  std::vector<int32_t> parent;
  std::vector<int32_t> child;
  int32_t root;
  bool rooted;
};

inline NNITree nni_tree(const IntegerMatrix& edge) {
  const int32_t n_edge = int32_t(edge.nrow());
  NNITree tree;
  tree.parent.resize(n_edge);
  tree.child.resize(n_edge);
  for (int32_t i = 0; i < n_edge; ++i) {
    tree.parent[i] = static_cast<int32_t>(edge(i, 0));
    tree.child[i] = static_cast<int32_t>(edge(i, 1));
  }
  tree.root = tree.parent[n_edge - 1];
  tree.rooted = tree.parent[n_edge - 3] != tree.root;
  return tree;
}

// Node (0-based) whose split duplicates the root split in a rooted tree
inline int32_t nni_trivial_two(const NNITree& tree, const bool rooted) {
  return rooted ? tree.child.back() - 1 : NOT_TRIVIAL;
}

inline int32_t nni_n_splits(const int32_t n_edge, const int32_t n_tip,
                            const bool rooted) {
  return n_edge - (rooted ? 1 : 0) - n_tip;
}

inline void nni_edge_to_splits(const NNITree& tree,
                               const int32_t n_tip,
                               const int32_t n_bin,
                               const int32_t trivial_two,
                               uint64_t* splits,
                               int32_t* names) {

  ASSERT((n_bin == n_tip + SL_BIN_SIZE - 1) / SL_BIN_SIZE);
  const int32_t n_edge = int32_t(tree.parent.size());
  const int32_t n_node = n_edge + 1;
  const int32_t trivial_origin = tree.root - 1;
  std::vector<uint64_t> tmp_splits(n_node * n_bin, 0);

  for (int32_t i = 0; i < n_tip; ++i) {
    const size_t bin_idx = static_cast<size_t>(i / SL_BIN_SIZE);
    const size_t idx = i * n_bin + bin_idx;
    ASSERT(idx < tmp_splits.size());
    ASSERT((i % SL_BIN_SIZE) < 64); // 1 << 64 = 0
    tmp_splits[idx] = static_cast<uint64_t>(1) << (i % SL_BIN_SIZE);
  }

  for (int32_t i = 0; i < n_edge - 1; ++i) { /* final edge is second root edge */
    const ptrdiff_t parent_row_i = (tree.parent[i] - 1) * n_bin;
    const ptrdiff_t child_row_i = (tree.child[i] - 1) * n_bin;
    for (int32_t j = 0; j < n_bin; ++j) {
      tmp_splits[parent_row_i + j] |= tmp_splits[child_row_i + j];
    }
  }

  int32_t n_trivial = 0;
  for (int32_t i = n_tip; i < n_node; ++i) {
    if (i == trivial_origin || i == trivial_two) {
      ++n_trivial;
    } else {
      const ptrdiff_t i_row = i * n_bin;
      for (int32_t j = 0; j < n_bin; ++j) {
        const size_t idx = i - n_tip - n_trivial;
        ASSERT(idx >= 0 && idx < n_splits);
        splits[idx * n_bin + j] = tmp_splits[i_row + j];
        names[idx] = i + 1;
      }
    }
  }
}

void nni_rf_matching (
    const uint64_t* a,
    const uint64_t* b,
    const int32_t n_splits,
    const int32_t n_bins,
    const int32_t n_tips,
    grf_match& matching) {

    ASSERT(n_tips > 3);
    ASSERT(n_tips <= NNI_MAX_TIPS); // Validated by exported wrapper

    const int32_t last_bin = n_bins - 1;
    const int32_t unset_tips = (n_tips % SL_BIN_SIZE) ?
      SL_BIN_SIZE - n_tips % SL_BIN_SIZE : 0;

    const uint64_t unset_mask = ALL_ONES >> unset_tips;

    matching.assign(n_splits, NA_INT32);
    HybridBuffer<uint64_t, NNI_STACK_SPLITS * NNI_STACK_BINS>
      b_complement(n_splits * n_bins);

    for (int32_t i = 0; i < n_splits; ++i) {
      const int32_t row_i = i * n_bins;
      for (int32_t bin = 0; bin < last_bin; ++bin) {
        const int32_t cell = row_i + bin;
        b_complement[cell] = ~b[cell];
      }
      const int32_t last_cell = row_i + last_bin;
      b_complement[last_cell] = b[last_cell] ^ unset_mask;
    }

    for (int32_t ai = 0; ai < n_splits; ++ai) {
      const int32_t a_row = ai * n_bins;
      for (int32_t bi = 0; bi < n_splits; ++bi) {

        bool all_match = true;
        bool all_complement = true;
        const int32_t b_row = bi * n_bins;

        for (int32_t bin = 0; bin < n_bins; ++bin) {
          const int32_t a_cell = int32_t(a_row) + int32_t(bin);
          const int32_t b_cell = int32_t(b_row) + int32_t(bin);
          if (a[a_cell] != b[b_cell]) {
            all_match = false;
            break;
          }
        }
        if (!all_match) {
          for (int32_t bin = 0; bin < n_bins; ++bin) {
            const int32_t a_cell = int32_t(a_row) + int32_t(bin);
            const int32_t bc_cell = int32_t(b_row) + int32_t(bin);
            if (a[a_cell] != b_complement[bc_cell]) {
              all_complement = false;
              break;
            }
          }
        }
        if (all_match || all_complement) {
          matching[ai] = bi + 1;
          break; // Only one match possible per split
        }
      }
    }
  }

// Writes lower, best_upper, tight_upper, loose_upper, fack_upper and li_upper
// to `bounds`.  `names_1` gives the node that defines each split of tree1;
// `match` maps these splits to those of the other tree.
void nni_score(const NNITree& tree1, const int32_t* names_1,
               const grf_match& match, const int32_t n_tip,
               int* bounds) {
  const int32_t node_0 = n_tip;
  const int32_t node_0_r = n_tip + 1;
  const int32_t n_edge = int32_t(tree1.parent.size());
  const bool rooted = tree1.rooted;

  int32_t lower_bound = 0;
  int32_t tight_score_bound = 0;

  int32_t loose_score_bound = 0;
  int32_t best_upper_bound = 0;
  int32_t fack_score_bound = 0;
  int32_t li_score_bound = 0;

  const int32_t n_distinct_edge = int32_t(n_edge - (rooted ? 1 : 0));

  HybridBuffer<bool, NNI_STACK_TIPS> matched_1(n_tip);
  HybridBuffer<int32_t, NNI_STACK_TIPS> unmatched_below(n_tip);

  const int32_t match_size = static_cast<int32_t>(match.size());
  for (int32_t i = 0; i < match_size; ++i) {
    ASSERT(n_edge != n_tip && n_tip > 3); // else names_1 uninitialized
    int32_t node_i = names_1[i] - node_0_r;
    if (match[i] == NA_INT32) {
      matched_1[node_i] = false;
      unmatched_below[node_i] = 1;
    } else {
      matched_1[node_i] = true;
    }
  }

  const int32_t edges_to_check = n_distinct_edge - (rooted ? 1 : 0);
  for (int32_t i = 0; i < edges_to_check; ++i) {
    const int32_t parent_i = tree1.parent[i] - 1;
    const int32_t child_i = tree1.child[i] - 1;
    // If edge is unmatched, add one to subtree size.
    if (child_i >= n_tip) {
      if (!matched_1[child_i - node_0]) {
        unmatched_below[parent_i - node_0] += unmatched_below[child_i - node_0];
      } else {
        update_score(unmatched_below[child_i - node_0], &lower_bound,
                     &tight_score_bound, &best_upper_bound, &loose_score_bound,
                     &li_score_bound, &fack_score_bound);
      }
    }
  }

  // Root edges:
  const int32_t root_node = tree1.root - node_0_r;

  if (rooted) {
    const int32_t root_child_1 = tree1.child[n_edge - 1] - 1;
    const int32_t root_child_2 = tree1.child[n_edge - 2] - 1;
    const int32_t unmatched_1 = root_child_1 < n_tip ? 0 :
                       unmatched_below[root_child_1 - node_0];
    if (root_child_2 >= n_tip) {
      const int32_t unmatched_2 = (root_child_2 < n_tip ? 0 :
                                   unmatched_below[root_child_2 - node_0]);
      if (!matched_1[root_child_2 - node_0]) {
        update_score(int32_t(unmatched_below[root_node]
                             + unmatched_1
                             + unmatched_2), &lower_bound,
                     &tight_score_bound, &best_upper_bound, &loose_score_bound,
                     &li_score_bound, &fack_score_bound);
      } else {
        update_score(unmatched_1, &lower_bound,
                     &tight_score_bound, &best_upper_bound, &loose_score_bound,
                     &li_score_bound, &fack_score_bound);
        update_score(unmatched_2, &lower_bound,
                     &tight_score_bound, &best_upper_bound, &loose_score_bound,
                     &li_score_bound, &fack_score_bound);
      }
    } else {
      update_score(unmatched_1, &lower_bound,
                   &tight_score_bound, &best_upper_bound, &loose_score_bound,
                   &li_score_bound, &fack_score_bound);
    }
  } else {
    update_score(unmatched_below[root_node], &lower_bound,
                 &tight_score_bound, &best_upper_bound, &loose_score_bound,
                 &li_score_bound, &fack_score_bound);
  }

  bounds[0] = lower_bound;
  bounds[1] = best_upper_bound;
  bounds[2] = tight_score_bound == NA_INT32 ? NA_INTEGER : tight_score_bound;
  bounds[3] = loose_score_bound;
  bounds[4] = fack_score_bound;
  bounds[5] = li_score_bound;
}

inline CharacterVector nni_bound_names() {
  return CharacterVector::create("lower", "best_upper", "tight_upper",
                                 "loose_upper", "fack_upper", "li_upper");
}

// [[Rcpp::export]]
IntegerVector cpp_nni_distance(const IntegerMatrix& edge1,
                               const IntegerMatrix& edge2,
                               const IntegerVector& nTip) {

  ASSERT(nTip[0] <= NNI_MAX_TIPS && "Cannot calculate NNI distance for trees with so many tips.");
  const int32_t n_tip = static_cast<int32_t>(nTip[0]);
  const int32_t n_edge = int32_t(edge1.nrow());

  ASSERT(n_edge == int32_t(edge2.nrow())
         && "Both trees must have the same number of edges.");

  IntegerVector ret(N_NNI_BOUNDS);
  ret.names() = nni_bound_names();
  if (n_tip < 4) {
    return ret;
  }

  const NNITree tree1 = nni_tree(edge1);
  const NNITree tree2 = nni_tree(edge2);
  const bool rooted = tree1.rooted;
  const int32_t n_bin = static_cast<int32_t>(((n_tip - 1) / SL_BIN_SIZE) + 1);
  const int32_t n_splits = nni_n_splits(n_edge, n_tip, rooted);

  if (n_splits < 1) {
    ASSERT(false && "NNI distance is undefined for trees with no splits"); // #nocov
  }

  std::unique_ptr<uint64_t[]> splits1(new uint64_t[n_splits * n_bin]);
  std::unique_ptr<uint64_t[]> splits2(new uint64_t[n_splits * n_bin]);
  std::unique_ptr<int32_t[]> names_1(new int32_t[n_splits]);
  ensure_cache();

  if (n_edge != n_tip) {
    nni_edge_to_splits(tree2, n_tip, n_bin, nni_trivial_two(tree2, rooted),
                       splits2.get(), names_1.get());
    nni_edge_to_splits(tree1, n_tip, n_bin, nni_trivial_two(tree1, rooted),
                       splits1.get(), names_1.get());
  } // else no internal nodes resolved

  grf_match match;
  nni_rf_matching(splits1.get(), splits2.get(), n_splits, n_bin, n_tip, match);
  nni_score(tree1, names_1.get(), match, n_tip, INTEGER(ret));

  return ret;
}

// One column per pair of trees, in the order of a `dist` object.
// Edges must be listed in postorder, with leaves numbered consistently.
// [[Rcpp::export]]
IntegerMatrix cpp_nni_distance_all_pairs(const List& edges,
                                         const IntegerVector& nTip) {
  const int32_t n_tip = static_cast<int32_t>(nTip[0]);
  ASSERT(n_tip <= NNI_MAX_TIPS);
  const R_xlen_t n_tree = edges.size();
  const R_xlen_t n_pair = n_tree * (n_tree - 1) / 2;

  IntegerMatrix ret(N_NNI_BOUNDS, n_pair);
  rownames(ret) = nni_bound_names();
  if (n_tip < 4 || n_pair == 0) {
    return ret;
  }
  ensure_cache();

  std::vector<NNITree> trees;
  trees.reserve(n_tree);
  for (R_xlen_t i = 0; i < n_tree; ++i) {
    trees.push_back(nni_tree(as<IntegerMatrix>(edges[i])));
  }

  const int32_t n_edge = int32_t(trees[0].parent.size());
  const int32_t n_bin = static_cast<int32_t>(((n_tip - 1) / SL_BIN_SIZE) + 1);
  const size_t max_splits = size_t(std::max(n_edge - n_tip, 0));
  const size_t stride = max_splits * n_bin;

  // Splits of each tree, computed once under its own rootedness
  std::vector<uint64_t> splits(stride * n_tree);
  std::vector<int32_t> names(max_splits * n_tree);
  if (n_edge != n_tip) {
    for (R_xlen_t i = 0; i < n_tree; ++i) {
      nni_edge_to_splits(trees[i], n_tip, n_bin,
                         nni_trivial_two(trees[i], trees[i].rooted),
                         &splits[stride * i], &names[max_splits * i]);
    }
  }

  // A pair is scored under the first tree's rootedness, so a second tree
  // that differs needs its splits recomputing
  std::vector<uint64_t> recast(stride);
  std::vector<int32_t> recast_names(max_splits);

  grf_match match;
  int* out = INTEGER(ret);
  R_xlen_t pair = 0;
  for (R_xlen_t i = 0; i < n_tree - 1; ++i) {
    const NNITree& tree_i = trees[i];
    const int32_t n_splits = nni_n_splits(n_edge, n_tip, tree_i.rooted);
    const uint64_t* splits_i = &splits[stride * i];
    for (R_xlen_t j = i + 1; j < n_tree; ++j) {
      const uint64_t* splits_j = &splits[stride * j];
      if (trees[j].rooted != tree_i.rooted) {
        nni_edge_to_splits(trees[j], n_tip, n_bin,
                           nni_trivial_two(trees[j], tree_i.rooted),
                           recast.data(), recast_names.data());
        splits_j = recast.data();
      }
      nni_rf_matching(splits_i, splits_j, n_splits, n_bin, n_tip, match);
      nni_score(tree_i, &names[max_splits * i], match, n_tip,
                out + pair * N_NNI_BOUNDS);
      if (++pair % 4096 == 0) {
        Rcpp::checkUserInterrupt();
      }
    }
  }

  return ret;
}
