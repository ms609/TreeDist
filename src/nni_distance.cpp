#include <math.h>
#include <Rcpp.h>
#include "tree_distances.h"
#include "splits.h"
#include "SplitList.h"

using namespace Rcpp;

const unsigned int NNI_MAX_TIPS = 5000; /* To avoid variable length arrays */

/* Exact value of diameter for trees with 0..N_EXACT edges, 
 * calculated by Li et al. 1996. */
const unsigned int N_EXACT = 11;
const unsigned int exact_diameter[] = {0, 0, 0, 0, 1, 3, 5, 7, 10, 12, 15, 18};

unsigned int sorting_number (unsigned int n_tip) {
  int log_ceiling = ceil(log2(n_tip));
  return (n_tip * log_ceiling - pow(2, log_ceiling) + 1);
}

unsigned int degenerate_distance (unsigned int n_tip) {
  /* Calculate the shortest length of the longest path in a tree.
   * To make this path as short as possible, divide tips into three 
   * balanced trees, joined by a single node that will form part of every 
   * longest path.  One of these subtrees will be filled with >= n/3 nodes */
  unsigned int nodes_in_full = (ceil(log2(n_tip / 3)) > 0) ? ceil(log2(n_tip / 3)) : 0U;
  /* We want to put a power of two tips in this subtree, such that every node is 
   * equally close to its root */
  unsigned int tips_in_full = pow(2, nodes_in_full);
  /* Now the remaining tips must be spread sub-evenly between the remaining 
   * edges from this node.  Picture halving the tips; removing tips from one side
   * until it is a power of two will reduce the number of nodes by one, whilst 
   * at worst (if this brings the other side over a power of two) increasing the 
   * other side by one. */
  unsigned int tips_left = n_tip - tips_in_full;
  unsigned int min_backbone_nodes = (nodes_in_full +
                                     ceil(log2(tips_left / 2)) + 1 > 0) ?
                                     nodes_in_full +
                                     ceil(log2(tips_left / 2)) + 1 : 0U;
  /* The worst-case scenario requires a move for every node not on the backbone: */
  unsigned int n_node = (n_tip > 2) ? n_tip - 2 : 0U;
  return (n_node - min_backbone_nodes);
}

void update_score (unsigned int *subtree_edges, int *tight_bound, 
                   int *loose_bound) {
  // Score subtree, add to score, and reset subtree size
  if (*tight_bound != NA_INTEGER && *subtree_edges <= N_EXACT) {
    *tight_bound += exact_diameter[*subtree_edges];
  } else {
    *tight_bound = NA_INTEGER;
  }
  *loose_bound += sorting_number(*subtree_edges + 3) +
    (2 * degenerate_distance(*subtree_edges + 2));
  
  *subtree_edges = 0;
}

// [[Rcpp::export]]
List cpp_nni_distance (IntegerMatrix edge1, IntegerMatrix edge2,
                       IntegerVector nTip) {
  unsigned int n_tip = nTip[0],
               max_nodes = n_tip + n_tip - 1,
               node_0 = n_tip + 1,
               subtree_size = 0U;
  int lower_bound = 0, tight_score_bound_1 = 0, loose_score_bound_1 = 0,
      tight_score_bound_2 = 0, loose_score_bound_2 = 0;
  if (n_tip > NNI_MAX_TIPS) {
    throw std::length_error("Cannot calculate NNI distance for trees with so many tips.");
  }
  
  RawMatrix splits1 = cpp_edge_to_splits(edge1, nTip);
  RawMatrix splits2 = cpp_edge_to_splits(edge2, nTip);
  CharacterVector names1 = rownames(splits1),
                  names2 = rownames(splits2);

  List matches = cpp_robinson_foulds_distance(splits1, splits2, nTip);
  
  IntegerVector match = matches[1];
  bool matched1[NNI_MAX_TIPS] {0}, matched2[NNI_MAX_TIPS] {0};
  for (unsigned int i = 0U; i != max_nodes - node_0; i++) {
    matched2[i] = false;
  }
  for (int i = 0; i != match.size(); i++) {
    if (match[i] == NA_INTEGER) {
      matched1[atoi(names1[i]) - node_0] = false;
      /* TODO counting here and not for matched2 assumes that trees are binary
       * (or resolution of tree1 > res(tree2)) */
      lower_bound++;
      Rcout << "Unmatched: " << names1[i] << "\n";
    } else {
      Rcout << "Matched: " << names1[i] << " - " << names2[match[i]] << "\n";
      matched1[atoi(names1[i]) - node_0] = true;
      matched2[atoi(names2[match[i]]) - node_0] = true;
    }
  }
  
  for (unsigned int i = 0U; i != (unsigned int) edge1.nrow(); i++) {
    unsigned int child_i = edge1(i, 1);
    // If edge is unmatched, add one to subtree size.
    if (child_i > n_tip && !matched1[child_i - node_0]) {
      subtree_size++;
    } else {
      Rcout << " Adding score for subtree [1] of size " << subtree_size << "\n";
      update_score(&subtree_size, &tight_score_bound_1, &loose_score_bound_1);
    }
  }
  update_score(&subtree_size, &tight_score_bound_1, &loose_score_bound_1);
  Rcout << "Score bounds from 1: " << tight_score_bound_1 << ", " 
        << loose_score_bound_1 << ".\n";
  
  
  for (unsigned int i = 0U; i != (unsigned int) edge2.nrow(); i++) {
    int child_i = edge2(i, 1);
    // If edge is unmatched, add one to subtree size.
    if (child_i > n_tip && !matched2[child_i - node_0]) {
      subtree_size++;
    } else {
      Rcout << " Adding score for subtree [2] of size " << subtree_size << "\n";
      update_score(&subtree_size, &tight_score_bound_2, &loose_score_bound_2);
    }
  }
  update_score(&subtree_size, &tight_score_bound_2, &loose_score_bound_2);
  Rcout << "Score bounds from 2: " << tight_score_bound_2 << ", " 
        << loose_score_bound_2 << ".\n";
  
  List ret = List::create(Named("lower") = lower_bound, 
                                _["tight_upper"] = tight_score_bound_1,
                                _["loose_upper"] = loose_score_bound_1);
  return (ret);
}
