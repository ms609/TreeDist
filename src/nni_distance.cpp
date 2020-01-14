#include <cmath>
#include <Rcpp.h>
#include "tree_distances.h"
#include "splits.h"
#include "SplitList.h"

using namespace Rcpp;

const unsigned int NNI_MAX_TIPS = 5000; /* To avoid variable length arrays */

/* Exact value of diameter for trees with 0..N_EXACT edges, 
 * calculated by Li et al. 1996. */
const unsigned int N_EXACT = 12;
const unsigned int exact_diameter[] = {0, 0, 0, 0, 1, 3, 5, 7, 10, 12, 15, 18, 21};

unsigned int sorting_number (unsigned int n_tip) {
  int log_ceiling = std::ceil(log2(n_tip));
  return (n_tip * log_ceiling - std::pow(2, log_ceiling) + 1);
}

unsigned int degenerate_distance (unsigned int n_tip) {
  /* Calculate the shortest length of the longest path in a tree.
   * To make this path as short as possible, divide tips into three 
   * balanced trees, joined by a single node that will form part of every 
   * longest path.  One of these subtrees will be filled with >= n/3 nodes */
  unsigned int nodes_in_full = (std::ceil(log2(n_tip / 3)) > 0) ? 
  std::ceil(log2(n_tip / 3)) : 0U;
  /* We want to put a power of two tips in this subtree, such that every node is 
   * equally close to its root */
  unsigned int tips_in_full = std::pow(2, nodes_in_full);
  /* Now the remaining tips must be spread sub-evenly between the remaining 
   * edges from this node.  Picture halving the tips; removing tips from one side
   * until it is a power of two will reduce the number of nodes by one, whilst 
   * at worst (if this brings the other side over a power of two) increasing the 
   * other side by one. */
  unsigned int tips_left = n_tip - tips_in_full;
  unsigned int min_backbone_nodes = 
    (nodes_in_full + std::ceil(log2(tips_left / 2)) + 1 > 0) ?
     nodes_in_full + std::ceil(log2(tips_left / 2)) + 1 : 0U;
  /* The worst-case scenario requires a move for every node not on the backbone: */
  unsigned int n_node = (n_tip > 2) ? n_tip - 2 : 0U;
  return (n_node - min_backbone_nodes);
}

// Score subtree, add to score, and reset subtree size
void update_score (const unsigned int subtree_edges, int *tight_bound, 
                   int *loose_bound) {
  if (subtree_edges) {
/*    Rcout << "   Region of " << subtree_edges << " unmatched edges.\n";*/
    const unsigned int subtree_tips = subtree_edges + 3;
    if (*tight_bound != NA_INTEGER && subtree_tips <= N_EXACT) {
/*      Rcout << "     Tight (" << *tight_bound << ") += " << exact_diameter[subtree_tips] << "\n";*/
      *tight_bound += exact_diameter[subtree_tips];
    } else {
      *tight_bound = NA_INTEGER;
    }
    /*Rcout << "     Loose (" << *loose_bound << ") += sorting No[" << subtree_tips << "], " << sorting_number(subtree_tips) 
            << "\n              +  2xO(" << subtree_tips << ") = 2x" << degenerate_distance(subtree_tips) << "\n";*/
    *loose_bound += sorting_number(subtree_tips) +
      (2 * degenerate_distance(subtree_tips));
  }
}

// [[Rcpp::export]]
IntegerVector cpp_nni_distance (IntegerMatrix edge1, IntegerMatrix edge2,
                       IntegerVector nTip) {
  unsigned int n_tip = nTip[0],
               node_0 = n_tip,
               node_0_r = n_tip + 1,
               n_edge = edge1.nrow();
  int lower_bound = 0, tight_score_bound = 0, loose_score_bound = 0;
  if (n_tip < 4) {
    return(IntegerVector::create(Named("lower") = 0,
                                 _["tight_upper"] = 0,
                                 _["loose_upper"] = 0));
  }
  if (n_tip > NNI_MAX_TIPS) {
    throw std::length_error("Cannot calculate NNI distance for trees with "
                              "so many tips.");
  }
  
  RawMatrix splits1 = cpp_edge_to_splits(edge1, nTip);
  RawMatrix splits2 = cpp_edge_to_splits(edge2, nTip);
  CharacterVector names1 = rownames(splits1);
  
  std::vector<int> match = cpp_robinson_foulds_distance(splits1, splits2, nTip)[1];
  
  bool matched1[NNI_MAX_TIPS] = {0};
  unsigned int unmatched_below[NNI_MAX_TIPS] = {0};

  for (int i = 0; i != int(match.size()); i++) {
    unsigned int node_i = atoi(names1[i]) - node_0_r;
    if (match[i] == NA_INTEGER) {
      matched1[node_i] = false;
      unmatched_below[node_i] = 1;
      lower_bound++;
    } else {
      matched1[node_i] = true;
    }
  }
  
  for (unsigned int i = 0U; i != (unsigned int) n_edge - 2; i++) {
    const unsigned int parent_i = edge1(i, 0) - 1, child_i = edge1(i, 1) - 1;
    // If edge is unmatched, add one to subtree size.
    if (child_i >= n_tip) {
      if (!matched1[child_i - node_0]) {
/*        Rcout << " Unmatched: " << (parent_i+1) << "-" << (child_i+1)<<"\n";*/
        unmatched_below[parent_i - node_0] += unmatched_below[child_i - node_0];
      } else {
/*        Rcout << " Matched: " << (parent_i+1) << "-" << (child_i+1) <<"\n";*/
        update_score(unmatched_below[child_i - node_0],
                     &tight_score_bound, &loose_score_bound);
      }
    }
  }
  // Root edges:
  const unsigned int root_node = edge1(n_edge - 2, 0) - node_0_r,
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

  return IntegerVector::create(Named("lower") = lower_bound, 
                               _["tight_upper"] = tight_score_bound,
                               _["loose_upper"] = loose_score_bound);
  
}
