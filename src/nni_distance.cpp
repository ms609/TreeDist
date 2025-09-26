#include <Rcpp/Lightest>
#include <TreeTools/SplitList.h>
#include <TreeTools/assert.h>

#include <algorithm> // for fill_n
#include <cmath>
#include <memory> // for unique_ptr
#include <mutex> // for once_flag

#include "tree_distances.h"
#include "li_diameters.h"

using namespace Rcpp;
using TreeTools::SplitList;

constexpr int32_t NA_INT32 = std::numeric_limits<int32_t>::min();


template <typename T, std::size_t StackSize>
class HybridBuffer {
public:
  explicit HybridBuffer(std::size_t n) : n_(n), data_(nullptr) {
    if (n_ <= StackSize) {
      data_ = stack_;
    } else {
      heap_ = std::unique_ptr<T[]>(new T[n_]());                       // #nocov
      data_ = heap_.get();                                             // #nocov
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
  alignas(T) T stack_[StackSize]{};
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
  
  for (int32_t i = 4; i != NNI_MAX_TIPS + 1; i++) {
    fack_lookup[i] = int32_t(((i - 2 - 2) * lg2_ceiling[i - 2]) + i - 2);
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
                  int32_t *closest_lower_bound,
                  int32_t *tight_bound, 
                  int32_t *closest_upper_bound,
                  int32_t *loose_bound, 
                  int32_t *li_bound, int32_t *fack_bound) {
  if (subtree_edges) {
    const int32_t subtree_tips = subtree_edges + 3;
    const int32_t lower = min_diameter[subtree_edges];
    const int32_t this_li = li[subtree_tips],
      fack = fack_lookup[subtree_tips],
      upper = this_li < fack ? this_li : fack;
    
    *lower_bound += lower;
    *li_bound += this_li;
    *fack_bound += fack;
    *loose_bound += upper;
    
    if (subtree_tips <= N_EXACT) {
      const int32_t exact = exact_diameter[subtree_tips];
      *closest_lower_bound += exact;
      if (*tight_bound != NA_INT32) {
        *tight_bound += exact;
      }
      *closest_upper_bound += exact;
    } else {
      *closest_lower_bound += lower;
      *tight_bound = NA_INT32;
      *closest_upper_bound += upper;
    }
  }
}


// Edges must be listed in postorder
inline void nni_edge_to_splits(const IntegerMatrix& edge,
                               const int32_t n_tip,
                               const int32_t n_edge,
                               const int32_t n_node,
                               const int32_t n_bin,
                               const int32_t trivial_origin,
                               const int32_t trivial_two,
                               std::unique_ptr<uint64_t[]>& splits,
                               std::unique_ptr<int32_t[]>& names) {
  
  ASSERT((n_bin == n_tip + SL_BIN_SIZE - 1) / SL_BIN_SIZE);
  std::vector<uint64_t> tmp_splits(n_node * n_bin, 0);
  
  for (int32_t i = 0; i < n_tip; ++i) {
    const size_t bin_idx = static_cast<size_t>(i / SL_BIN_SIZE);
    const size_t idx = i * n_bin + bin_idx;
    ASSERT(idx < tmp_splits.size());
    ASSERT((i % SL_BIN_SIZE) < 64); // 1 << 64 = 0
    tmp_splits[idx] = static_cast<uint64_t>(1) << (i % SL_BIN_SIZE);
  }
  
  for (int32_t i = 0; i < n_edge - 1; ++i) { /* final edge is second root edge */
    const ptrdiff_t parent_row_i = (edge(i, 0) - 1) * n_bin;
    const ptrdiff_t child_row_i = (edge(i, 1) - 1) * n_bin;
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
  
grf_match nni_rf_matching (
    const std::unique_ptr<uint64_t[]>& a, 
    const std::unique_ptr<uint64_t[]>& b,
    const int32_t n_splits,
    const int32_t n_bins,
    const int32_t n_tips) {
  
    ASSERT(n_splits > 0);
    ASSERT(n_tips > 3);
    if (n_tips > NNI_MAX_TIPS) {
      Rcpp::stop("Cannot calculate NNI distance for trees with so many tips.");
    }
    
    // #nocov start
    if (static_cast<int64_t>(n_splits) * static_cast<int64_t>(n_bins) >
          static_cast<int64_t>(std::numeric_limits<int32_t>::max())) {
      Rcpp::stop("Cannot calculate NNI distance for trees with so many splits.");
    }
    // #nocov end
    
    const int32_t last_bin = n_bins - 1;
    const int32_t unset_tips = (n_tips % SL_BIN_SIZE) ? 
      SL_BIN_SIZE - n_tips % SL_BIN_SIZE : 0;
    
    const uint64_t unset_mask = ALL_ONES >> unset_tips;
    
    grf_match matching(n_splits, NA_INT32);
    for (int32_t i = 0; i != n_splits; i++) {
      ASSERT(matching[i] == NA_INT32);
    }
    
    HybridBuffer<uint64_t, NNI_STACK_SPLITS * NNI_STACK_BINS>
      b_complement(n_splits * n_bins);
    
    for (int32_t i = 0; i < n_splits; ++i) {
      const int32_t row_i = i * n_bins;
      for (int32_t bin = 0; bin < last_bin; ++bin) {
        const int32_t cell = row_i + bin;
        ASSERT(cell >= 0 && cell < int32_t(n_splits) * int32_t(n_bins));
        b_complement[cell] = ~b[cell];
      }
      const int32_t last_cell = row_i + last_bin;
      ASSERT(last_cell >= 0 && last_cell < int32_t(n_splits) * int32_t(n_bins));
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
          ASSERT(a_cell >= 0 && a_cell < int32_t(n_splits) * int32_t(n_bins));
          ASSERT(b_cell >= 0 && b_cell < int32_t(n_splits) * int32_t(n_bins));
          if (a[a_cell] != b[b_cell]) {
            all_match = false;
            break;
          }
        }
        if (!all_match) {
          for (int32_t bin = 0; bin < n_bins; ++bin) {
            const int32_t a_cell = int32_t(a_row) + int32_t(bin);
            const int32_t bc_cell = int32_t(b_row) + int32_t(bin);
            ASSERT(a_cell >= 0 && a_cell < int32_t(n_splits) * int32_t(n_bins));
            ASSERT(bc_cell >= 0 && bc_cell < int32_t(n_splits) * int32_t(n_bins));
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
    
    return matching;
  }

// [[Rcpp::export]]
IntegerVector cpp_nni_distance(const IntegerMatrix& edge1, 
                               const IntegerMatrix& edge2,
                               const IntegerVector& nTip) {
  
  if (nTip[0] > NNI_MAX_TIPS) {
    Rcpp::stop("Cannot calculate NNI distance for trees with "
                            "so many tips.");
  }
  const int32_t n_tip = static_cast<int32_t>(nTip[0]);
  const int32_t node_0 = n_tip;
  const int32_t node_0_r = n_tip + 1;
  const int32_t n_edge = int32_t(edge1.nrow());
  
  int32_t lower_bound = 0;
  int32_t best_lower_bound = 0;
  int32_t tight_score_bound = 0;
  
  int32_t loose_score_bound = 0;
  int32_t best_upper_bound = 0;
  int32_t fack_score_bound = 0;
  int32_t li_score_bound = 0;
  
  if (n_edge != int32_t(edge2.nrow())) {
    Rcpp::stop("Both trees must have the same number of edges. "
                            "Is one rooted and the other unrooted?");
  }

  if (n_tip < 4) {
    return(IntegerVector::create(Named("lower") = 0,
                                 _["best_lower"] = 0,
                                 _["tight_upper"] = 0,
                                 _["best_upper"] = 0,
                                 _["loose_upper"] = 0,
                                 _["fack_upper"] = 0,
                                 _["li_upper"] = 0));
  }
  
  const int32_t root_1 = static_cast<int32_t>(edge1(n_edge - 1, 0));
  const int32_t root_2 = static_cast<int32_t>(edge2(n_edge - 1, 0));
  bool rooted = static_cast<int32_t>(edge1(n_edge - 3, 0)) != root_1;
  constexpr int32_t NOT_TRIVIAL = std::numeric_limits<int32_t>::max();
  const int32_t n_node = n_edge + 1;
  const int32_t n_bin = static_cast<int32_t>(((n_tip - 1) / SL_BIN_SIZE) + 1);
  const int32_t trivial_origin_1 = root_1 - 1;
  const int32_t trivial_origin_2 = root_2 - 1;
  const int32_t trivial_two_1 = rooted ?
    static_cast<int32_t>(edge1(n_edge - 1, 1)) - 1 :
    NOT_TRIVIAL;
  const int32_t trivial_two_2 = rooted ?
    static_cast<int32_t>(edge2(n_edge - 1, 1)) - 1 :
    NOT_TRIVIAL;
  const int32_t n_distinct_edge = int32_t(n_edge - (rooted ? 1 : 0));
  const int32_t n_splits = n_distinct_edge - n_tip;
  
  if (n_splits < 1) {
    Rcpp::stop("NNI distance is undefined for trees with no splits"); // #nocov
  }
  
  std::unique_ptr<uint64_t[]> splits1(new uint64_t[n_splits * n_bin]);
  std::unique_ptr<uint64_t[]> splits2(new uint64_t[n_splits * n_bin]);
  // std::vector<int32_t> names_1;
  // names_1.reserve(n_splits);
  std::unique_ptr<int32_t[]> names_1(new int32_t[n_splits]);
  ensure_cache();
  
  if (n_edge != n_tip && n_tip > 3) {
    nni_edge_to_splits(edge2, n_tip, n_edge, n_node, n_bin, 
                       trivial_origin_2, trivial_two_2, splits2, names_1);
    nni_edge_to_splits(edge1, n_tip, n_edge, n_node, n_bin,
                       trivial_origin_1, trivial_two_1, splits1, names_1);
  } // else no internal nodes resolved
  
  grf_match match = nni_rf_matching(splits1, splits2, n_splits, n_bin, n_tip);
  
  HybridBuffer<bool, NNI_STACK_TIPS> matched_1(n_tip);
  HybridBuffer<int32_t, NNI_STACK_TIPS> unmatched_below(n_tip);
  
  const int32_t match_size = static_cast<int32_t>(match.size());
  for (int32_t i = 0; i < match_size; ++i) {
    ASSERT(n_edge != n_tip && n_tip > 3); // else names_1 uninitialized
    int32_t node_i = names_1[i] - node_0_r;
    ASSERT(node_i >= 0 && node_i < n_tip);
    if (match[i] == NA_INT32) {
      matched_1[node_i] = false;
      unmatched_below[node_i] = 1;
    } else {
      matched_1[node_i] = true;
    }
  }
  
  const int32_t edges_to_check = n_distinct_edge - (rooted ? 1 : 0);
  for (int32_t i = 0; i < edges_to_check; ++i) {
    const int32_t parent_i = static_cast<int32_t>(edge1(i, 0)) - 1;
    const int32_t child_i = static_cast<int32_t>(edge1(i, 1)) - 1;
    // If edge is unmatched, add one to subtree size.
    if (child_i >= n_tip) {
      if (!matched_1[child_i - node_0]) {
        unmatched_below[parent_i - node_0] += unmatched_below[child_i - node_0];
      } else {
        update_score(unmatched_below[child_i - node_0], &lower_bound,
                     &best_lower_bound, &tight_score_bound,
                     &best_upper_bound, &loose_score_bound,
                     &li_score_bound, &fack_score_bound);
      }
    }
  }
  
  // Root edges:
  const int32_t root_node = root_1 - node_0_r;
  
  if (rooted) {
    const int32_t root_child_1 = static_cast<int32_t>(edge1(n_edge - 1, 1)) - 1;
    const int32_t root_child_2 = static_cast<int32_t>(edge1(n_edge - 2, 1)) - 1;
    const int32_t unmatched_1 = root_child_1 < n_tip ? 0 :
                       unmatched_below[root_child_1 - node_0];
    if (root_child_2 >= n_tip) {
      const int32_t unmatched_2 = (root_child_2 < n_tip ? 0 :
                                   unmatched_below[root_child_2 - node_0]);
      if (!matched_1[root_child_2 - node_0]) {
        update_score(int32_t(unmatched_below[root_node]
                             + unmatched_1
                             + unmatched_2), &lower_bound, &best_lower_bound,
                     &tight_score_bound, &best_upper_bound, &loose_score_bound,
                     &li_score_bound, &fack_score_bound);
      } else {
        update_score(unmatched_1, &lower_bound, &best_lower_bound, 
                     &tight_score_bound, &best_upper_bound, &loose_score_bound,
                     &li_score_bound, &fack_score_bound);
        update_score(unmatched_2, &lower_bound, &best_lower_bound, 
                     &tight_score_bound, &best_upper_bound, &loose_score_bound,
                     &li_score_bound, &fack_score_bound);
      }
    } else {
      update_score(unmatched_1, &lower_bound, &best_lower_bound, 
                   &tight_score_bound, &best_upper_bound, &loose_score_bound, 
                   &li_score_bound, &fack_score_bound);
    }
  } else {
    update_score(unmatched_below[root_node], &lower_bound, &best_lower_bound,
                 &tight_score_bound, &best_upper_bound, &loose_score_bound,
                 &li_score_bound, &fack_score_bound);
  }

  return IntegerVector::create(
    Named("lower")   = lower_bound,
    _["best_lower"] = best_lower_bound,
    _["tight_upper"] = tight_score_bound == NA_INT32 
                     ? NA_INTEGER : tight_score_bound,
    _["best_upper"] = best_upper_bound,
    _["loose_upper"] = loose_score_bound,
    _["fack_upper"] = fack_score_bound,
    _["li_upper"] = li_score_bound);
}
