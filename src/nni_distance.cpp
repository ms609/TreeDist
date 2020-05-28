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
                        const uint16* n_tip,
                        const uint16* n_edge,
                        const uint16* n_node,
                        const uint16* n_bin,
                        const uint16* trivial_origin,
                        const uint16* trivial_two,
                        splitbit* splits,
                        uint16* names) {
  
  
  splitbit** tmp_splits = new splitbit*[*n_node];
  for (uint16 i = 0; i != *n_node; i++) {
    tmp_splits[i] = new splitbit[*n_bin]();
  }
  
  for (uint16 i = 0; i != *n_tip; i++) {
    tmp_splits[i][uint16(i / BIN_SIZE)] = powers_of_two[i % BIN_SIZE];
  }
  
  for (uint16 i = 0; i != *n_edge - 1; i++) { /* final edge is second root edge */
    for (uint16 j = 0; j != *n_bin; j++) {
      tmp_splits[uint16(edge(i, 0) - 1)][j] |= tmp_splits[uint16(edge(i, 1) - 1)][j];
    }
  }
  
  for (uint16 i = 0; i != *n_tip; i++) {
    delete[] tmp_splits[i];
  }
  
  uint16 n_trivial = 0;
  for (uint16 i = *n_tip; i != *n_node; i++) {
    if (i == *trivial_origin || i == *trivial_two) {
      n_trivial++;
    } else {
      for (uint16 j = 0; j != *n_bin; j++) {
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
    const uint16* n_splits,
    const uint16* n_bins,
    const uint16* n_tips) {
    
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
  const uint16 n_tip = nTip[0],
               node_0 = n_tip,
               node_0_r = n_tip + 1,
               n_edge = edge1.nrow();
  int16 lower_bound = 0, tight_score_bound = 0, loose_score_bound = 0;
  if (n_edge != edge2.nrow()) {
    throw std::length_error("Both trees must have the same number of edges. "
                            "Is one rooted and the other unrooted?");
  }

  if (n_tip < 4) {
    return(IntegerVector::create(Named("lower") = 0,
                                 _["tight_upper"] = 0,
                                 _["loose_upper"] = 0));
  }
  
  const uint16
    NOT_TRIVIAL = UINTX_MAX,
    
    n_node = n_edge + 1,
    n_bin = ((n_tip - 1) / BIN_SIZE) + 1,
    
    trivial_origin_1 = edge1(n_edge - 1, 0) - 1,
    trivial_origin_2 = edge2(n_edge - 1, 0) - 1,
    
    trivial_two_1 = (edge1(n_edge - 1, 0) == edge1(n_edge - 3, 0) ?
                     NOT_TRIVIAL : (edge1(n_edge - 1, 1) - 1L)),
    trivial_two_2 = (edge2(n_edge - 1, 0) == edge2(n_edge - 3, 0) ?
                     NOT_TRIVIAL : (edge2(n_edge - 1, 1) - 1L)),
    
    n_splits = n_edge - n_tip - (trivial_two_1 != NOT_TRIVIAL ? 1 : 0)
  ;
  
  splitbit 
    *splits1 = new splitbit[n_splits * n_bin],
    *splits2 = new splitbit[n_splits * n_bin];
  uint16 *names1 = new uint16[n_splits];
  
  if (n_edge != n_tip && n_tip > 3) {
    nni_edge_to_splits(edge2, &n_tip, &n_edge, &n_node, &n_bin, 
                       &trivial_origin_2, &trivial_two_2, splits2, names1);
    nni_edge_to_splits(edge1, &n_tip, &n_edge, &n_node, &n_bin,
                       &trivial_origin_1, &trivial_two_1, splits1, names1);
  } // else no internal nodes resolved
  
  rf_match match = nni_rf_matching(splits1, splits2, &n_splits, &n_bin, &n_tip);
  
  delete[] splits1;
  delete[] splits2;
  
  bool matched1[NNI_MAX_TIPS] = {0};
  int16 unmatched_below[NNI_MAX_TIPS] = {0};

  for (int16 i = 0; i != int16(match.size()); i++) {
    int16 node_i = names1[i] - node_0_r;
    if (match[i] == NA_INT16) {
      matched1[node_i] = false;
      unmatched_below[node_i] = 1;
      lower_bound++;
    } else {
      matched1[node_i] = true;
    }
  }
  delete[] names1;
  
  for (int16 i = 0; i != n_edge - 2; i++) {
    const int16 parent_i = edge1(i, 0) - 1, child_i = edge1(i, 1) - 1;
    // If edge is unmatched, add one to subtree size.
    if (child_i >= n_tip) {
      if (!matched1[child_i - node_0]) {
        unmatched_below[parent_i - node_0] += unmatched_below[child_i - node_0];
      } else {
        update_score(unmatched_below[child_i - node_0],
                     &tight_score_bound, &loose_score_bound);
      }
    }
  }
  // Root edges:
  const int16 root_node = edge1(n_edge - 2, 0) - node_0_r,
              root_child_1 = edge1(n_edge - 2, 1) - 1,
              root_child_2 = edge1(n_edge - 1, 1) - 1,
              unmatched_2 = (root_child_2 < n_tip ? 0 :
                               unmatched_below[root_child_2 - node_0]);
  if (root_child_1 >= n_tip) {
    if (!matched1[root_child_1 - node_0]) {
/*      Rcout << " Root edge unmatched\n";*/
      update_score(unmatched_below[root_node]
                     + unmatched_below[root_child_1 - node_0] + unmatched_2,
                   &tight_score_bound, &loose_score_bound);
    } else {
/*      Rcout << " Root edge matched\n";*/
      update_score(unmatched_below[root_child_1 - node_0],
                   &tight_score_bound, &loose_score_bound);
      update_score(unmatched_2,
                   &tight_score_bound, &loose_score_bound);
    }
  } else {
    update_score(unmatched_2,
                 &tight_score_bound, &loose_score_bound);
  }

  return IntegerVector::create(Named("lower")   = lower_bound,
                               _["tight_upper"] = tight_score_bound == NA_INT16 
                                 ? NA_INTEGER : tight_score_bound,
                               _["loose_upper"] = loose_score_bound);
  
}
