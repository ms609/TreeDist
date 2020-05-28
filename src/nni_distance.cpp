#include <cmath>
#include <Rcpp.h>
#include "tree_distances.h"
#include "SplitList.h"

using namespace Rcpp;

const int16 NNI_MAX_TIPS = 5000; /* To avoid variable length arrays */

/* Exact value of diameter for trees with 0..N_EXACT edges, 
 * calculated by Li et al. 1996. */
const int16 N_EXACT = 12;
const int16 exact_diameter[] = {0, 0, 0, 0, 1, 3, 5, 7, 10, 12, 15, 18, 21};

int16 sorting_number (int16 n_tip) {
  const int16 log_ceiling = std::ceil(log2(n_tip));
  return (n_tip * log_ceiling - std::pow(2, log_ceiling) + 1);
}

int16 degenerate_distance (const int16 n_tip) {
  /* Calculate the shortest length of the longest path in a tree.
   * To make this path as short as possible, divide tips into three 
   * balanced trees, joined by a single node that will form part of every 
   * longest path.  One of these subtrees will be filled with >= n/3 nodes */
  const int16 nodes_in_full = (std::ceil(log2(n_tip / 3)) > 0) ? 
                               std::ceil(log2(n_tip / 3)) : 0;
  /* We want to put a power of two tips in this subtree, such that every node is 
   * equally close to its root */
  const int16 tips_in_full = std::pow(2, nodes_in_full);
  /* Now the remaining tips must be spread sub-evenly between the remaining 
   * edges from this node.  Picture halving the tips; removing tips from one side
   * until it is a power of two will reduce the number of nodes by one, whilst 
   * at worst (if this brings the other side over a power of two) increasing the 
   * other side by one. */
  const int16 tips_left = n_tip - tips_in_full;
  const int16 min_backbone_nodes = 
    (nodes_in_full + std::ceil(log2(tips_left / 2)) + 1 > 0) ?
     nodes_in_full + std::ceil(log2(tips_left / 2)) + 1 : 0;

  /* The worst-case scenario requires a move for every node not on the backbone: */
  const int16 n_node = (n_tip > 2) ? n_tip - 2 : 0;
  return (n_node - min_backbone_nodes);
}

// Score subtree, add to score, and reset subtree size
void update_score (const int16 subtree_edges, int16 *tight_bound, 
                   int16 *loose_bound) {
  if (subtree_edges) {
    const int16 subtree_tips = subtree_edges + 3;
    if (*tight_bound != NA_INT16 && subtree_tips <= N_EXACT) {
      *tight_bound += exact_diameter[subtree_tips];
    } else {
      *tight_bound = NA_INT16;
    }

    *loose_bound += sorting_number(subtree_tips) +
      (2 * degenerate_distance(subtree_tips));
  }
}


// Edges must be listed in postorder
void nni_edge_to_splits(const IntegerMatrix edge,
                        const int16* n_tip,
                        const int16* n_edge,
                        const int16* n_node,
                        const int16* n_bin,
                        const int16* trivial_origin,
                        const int16* trivial_two,
                        splitbit* splits,
                        int16* names) {
  
  
  splitbit** tmp_splits = new splitbit*[*n_node];
  for (int16 i = 0; i != *n_node; i++) {
    tmp_splits[i] = new splitbit[*n_bin]();
  }
  
  for (int16 i = 0; i != *n_tip; i++) {
    tmp_splits[i][int16(i / BIN_SIZE)] = powers_of_two[i % BIN_SIZE];
  }
  
  for (int16 i = 0; i != *n_edge - 1; i++) { /* final edge is second root edge */
    for (int16 j = 0; j != *n_bin; j++) {
      tmp_splits[int16(edge(i, 0) - 1)][j] |= tmp_splits[int16(edge(i, 1) - 1)][j];
    }
  }
  
  for (int16 i = 0; i != *n_tip; i++) {
    delete[] tmp_splits[i];
  }
  
  int16 n_trivial = 0;
  for (int16 i = *n_tip; i != *n_node; i++) {
    if (i == *trivial_origin || i == *trivial_two) {
      n_trivial++;
    } else {
      for (int16 j = 0; j != *n_bin; j++) {
        splits[((i - *n_tip - n_trivial) * *n_bin) + j] = tmp_splits[i][j];
        names[i - *n_tip - n_trivial] = (i + 1);
      }
    }
    delete[] tmp_splits[i];
  }
  
  delete[] tmp_splits;
}

rf_match nni_rf_matching (
    const splitbit* a, 
    const splitbit* b,
    const int16* n_splits,
    const int16* n_bins,
    const int16* n_tips) {
    
    const int16
      last_bin = *n_bins - 1,
      unset_tips = (*n_tips % BIN_SIZE) ? BIN_SIZE - *n_tips % BIN_SIZE : 0;
    const splitbit unset_mask = ALL_ONES >> unset_tips;
    
    rf_match matching (*n_splits);
    for (int16 i = 0; i != *n_splits; i++) matching[i] = NA_INT16;
    
    splitbit b_complement[MAX_SPLITS][MAX_BINS];
    for (int16 i = 0; i != *n_splits; i++) {
      for (int16 bin = 0; bin != last_bin; bin++) {
        b_complement[i][bin] = ~b[i * *n_bins + bin];
      }
      b_complement[i][last_bin] = b[i * *n_bins + last_bin] ^ unset_mask;
    }
    
    for (int16 ai = 0; ai != *n_splits; ai++) {
      for (int16 bi = 0; bi != *n_splits; bi++) {
        
        bool all_match = true, all_complement = true;
        
        for (int16 bin = 0; bin != *n_bins; bin++) {
          if ((a[ai * *n_bins + bin] != b[bi * *n_bins + bin])) {
            all_match = false;
            break;
          }
        }
        if (!all_match) {
          for (int16 bin = 0; bin != *n_bins; bin++) {
            if ((a[ai * *n_bins + bin] != b_complement[bi][bin])) {
              all_complement = false;
              break;
            }
          }
        }
        if (all_match || all_complement) {
          matching[ai] = bi + 1;
          break; /* Only one match possible per split */
        }
      }
    }
    
    return (matching);
  }

// [[Rcpp::export]]
IntegerVector cpp_nni_distance (const IntegerMatrix edge1, 
                                const IntegerMatrix edge2,
                                const IntegerVector nTip) {
  
  if (nTip[0] > NNI_MAX_TIPS) {
    throw std::length_error("Cannot calculate NNI distance for trees with "
                            "so many tips.");
  }
  const int16 
    n_tip = nTip[0],
    node_0 = n_tip,
    node_0_r = n_tip + 1,
    n_edge = edge1.nrow()
  ;
  int16 lower_bound = 0, tight_score_bound = 0, loose_score_bound = 0;
  if (n_edge != int16(edge2.nrow())) {
    throw std::length_error("Both trees must have the same number of edges. "
                            "Is one rooted and the other unrooted?");
  }

  if (n_tip < 4) {
    return(IntegerVector::create(Named("lower") = 0,
                                 _["tight_upper"] = 0,
                                 _["loose_upper"] = 0));
  }
  
  const int16 
    root_1 = edge1(n_edge - 1, 0),
    root_2 = edge2(n_edge - 1, 0)
  ;
  
  bool rooted = edge1(n_edge - 3, 0) != root_1;
  
  const uint16
    NOT_TRIVIAL = UINTX_MAX;
  
  const int16
    n_node = n_edge + 1,
    n_bin = ((n_tip - 1) / BIN_SIZE) + 1,
    
    trivial_origin_1 = root_1 - 1,
    trivial_origin_2 = root_2 - 1,
    
    trivial_two_1 = (rooted ? (edge1(n_edge - 1, 1) - 1) : NOT_TRIVIAL),
    trivial_two_2 = (rooted ? (edge2(n_edge - 1, 1) - 1) : NOT_TRIVIAL),
    
    n_distinct_edge = n_edge - (rooted ? 1 : 0),
    n_splits = n_distinct_edge - n_tip
  ;
  
  splitbit 
    *splits1 = new splitbit[n_splits * n_bin],
    *splits2 = new splitbit[n_splits * n_bin];
  int16 *names_1 = new int16[n_splits];
  
  if (n_edge != n_tip && n_tip > 3) {
    nni_edge_to_splits(edge2, &n_tip, &n_edge, &n_node, &n_bin, 
                       &trivial_origin_2, &trivial_two_2, splits2, names_1);
    nni_edge_to_splits(edge1, &n_tip, &n_edge, &n_node, &n_bin,
                       &trivial_origin_1, &trivial_two_1, splits1, names_1);
  } // else no internal nodes resolved
  
  rf_match match = nni_rf_matching(splits1, splits2, &n_splits, &n_bin, &n_tip);
  
  delete[] splits1;
  delete[] splits2;
  
  bool matched_1[NNI_MAX_TIPS] = {0};
  int16 unmatched_below[NNI_MAX_TIPS] = {0};

  for (int16 i = 0; i != int16(match.size()); i++) {
    int16 node_i = names_1[i] - node_0_r;
    if (match[i] == NA_INT16) {
      matched_1[node_i] = false;
      unmatched_below[node_i] = 1;
      lower_bound++;
    } else {
      matched_1[node_i] = true;
    }
  }
  delete[] names_1;
  
  for (int16 i = 0; i != n_distinct_edge - (rooted ? 1 : 0); i++) {
    const int16 parent_i = edge1(i, 0) - 1, child_i = edge1(i, 1) - 1;
    // If edge is unmatched, add one to subtree size.
    if (child_i >= n_tip) {
      if (!matched_1[child_i - node_0]) {
        unmatched_below[parent_i - node_0] += unmatched_below[child_i - node_0];
      } else {
        update_score(unmatched_below[child_i - node_0],
                     &tight_score_bound, &loose_score_bound);
      }
    }
  }
  
  // Root edges:
  const int16 root_node = root_1 - node_0_r;
  
  if (rooted) {
    const int16
      root_child_1 = edge1(n_edge - 1, 1) - 1,
      root_child_2 = edge1(n_edge - 2, 1) - 1,
      
      unmatched_1 = root_child_1 < n_tip ? 0 :
                       unmatched_below[root_child_1 - node_0]
    ;
    if (root_child_2 >= n_tip) {
      const int16 unmatched_2 = (root_child_2 < n_tip ? 0 :
                                   unmatched_below[root_child_2 - node_0]);
      if (!matched_1[root_child_2 - node_0]) {
        update_score(unmatched_below[root_node]
                       + unmatched_1
                       + unmatched_2,
                     &tight_score_bound, &loose_score_bound);
      } else {
        update_score(unmatched_1,
                     &tight_score_bound, &loose_score_bound);
        update_score(unmatched_2,
                     &tight_score_bound, &loose_score_bound);
      }
    } else {
      update_score(unmatched_1,
                   &tight_score_bound, &loose_score_bound);
    }
  } else {
    update_score(unmatched_below[root_node],
                   &tight_score_bound, &loose_score_bound);
  }

  return IntegerVector::create(Named("lower")   = lower_bound,
                               _["tight_upper"] = tight_score_bound == NA_INT16 
                                 ? NA_INTEGER : tight_score_bound,
                               _["loose_upper"] = loose_score_bound);
  
}
