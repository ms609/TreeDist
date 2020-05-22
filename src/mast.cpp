#include <cstdlib>
#include <Rcpp.h>
#include "ints.h"
using namespace Rcpp;

const int16 MAST_MAX_NODE = 700, /* Much larger --> overflow size of stack frame */
            MAST_MAX_TIP = MAST_MAX_NODE + 1,
            MAST_MAX_ALLNODE = MAST_MAX_NODE - 1 + MAST_MAX_NODE;

int16 max_of_six (const int16 m1,
                  const int16 m2,
                  const int16 m3,
                  const int16 m4,
                  const int16 m5,
                  const int16 m6) {
  int16 largest = (m1 > m2 ? m1 : m2);
  if (m3 > largest) largest = m3;
  if (m4 > largest) largest = m4;
  if (m5 > largest) largest = m5;
  return (m6 > largest ? m6 : largest);
}

void fill_M (const int16 node1, const int16 node2,
             int_fast16_t *M,
             const bool* t1_descendantsof,
             const bool* t2_descendantsof,
             const int16* t1_left,
             const int16* t1_right,
             const int16* t2_left,
             const int16* t2_right,
             const int16 n_tip,
             const int16 M_dim) {
  if (node1 < n_tip) { // Node 1 is a leaf
    if (node2 < n_tip) { // Node 1 and node 2 are leaves
      M[node1 * M_dim + node2] = node1 == node2 ? 1 : 0;
    } else { // Node 1 is a leaf; node 2 is internal
      M[node1 * M_dim + node2] = t2_descendantsof[(node2 - n_tip) * n_tip + node1] ? 1 : 0; 
    }
  } else {
    if (node2 < n_tip) { // Node 1 is internal; node 2 is a leaf
      M[node1 * M_dim + node2] = t1_descendantsof[(node1 - n_tip) * n_tip + node2] ? 1 : 0;
    } else { // Node 1 & Node 2 are internal
      const int16 l1 = t1_left[node1 - n_tip],
                  r1 = t1_right[node1 - n_tip],
                  l2 = t2_left[node2 - n_tip],
                  r2 = t2_right[node2 - n_tip];
      M[node1 * M_dim + node2] = max_of_six(
        M[l1 * M_dim + l2] + M[r1 * M_dim + r2],
        M[l1 * M_dim + r2] + M[r1 * M_dim + l2],
        M[node1 * M_dim + l2],
        M[node1 * M_dim + r2],
        M[l1 * M_dim + node2],
        M[r1 * M_dim + node2]);
    }
  }
}

// edge1 and edge2 should be edge matrices taken from trees of class `phylo`,
// with one subtracted from all entries so that leaves are numbered from 
// 0..(nTip - 1).
// [[Rcpp::export]]
int cpp_mast (IntegerMatrix edge1, IntegerMatrix edge2, IntegerVector nTip) {
  const int16 n_tip = nTip[0],
              n_internal = n_tip - 1,
              n_all_nodes = n_tip + n_internal,
              n_edge = edge1.nrow();
  if (edge2.nrow() != n_edge) {
    throw std::length_error("Both trees must contain the same number of edges.");
  }
  if (n_tip > MAST_MAX_TIP) {
    throw std::length_error("Tree too large to analyse on stack; please "
                              "contact maintainer for advice.");
  }
  
  int16
    t1_left[MAST_MAX_NODE] = {}, t1_right[MAST_MAX_NODE] = {},
    t2_left[MAST_MAX_NODE] = {}, t2_right[MAST_MAX_NODE] = {};
  bool t1_has_child[MAST_MAX_NODE] = {}, t2_has_child[MAST_MAX_NODE] = {};
  bool *t1_descendantsof = new bool[n_internal * n_tip](),
       *t2_descendantsof = new bool[n_internal * n_tip]();
  
  for (int16 i = 0; i != n_edge; i++) {
    const int16 
      parent1 = edge1(i, 0), child1 = edge1(i, 1), index1 = parent1 - n_tip,
      parent2 = edge2(i, 0), child2 = edge2(i, 1), index2 = parent2 - n_tip;
    if (t1_has_child[index1]) {
      t1_right[index1] = child1;
    } else {
      t1_left[index1] = child1;
      t1_has_child[index1] = true;
    }
    if (t2_has_child[index2]) {
      t2_right[index2] = child2; 
    } else {
      t2_left[index2] = child2;
      t2_has_child[index2] = true;
    }
    
    if (child1 < n_tip) { // Child 1 is a tip
      t1_descendantsof[index1 * n_tip + child1] = true;
    } else { // Child 1 is a node
      for (int16 tip = 0; tip != n_tip; tip++) {
        t1_descendantsof[index1 * n_tip + tip] |= 
          t1_descendantsof[(child1 - n_tip) * n_tip + tip];
      }
    }
    if (child2 < n_tip) {
      t2_descendantsof[index2 * n_tip + child2] = true;
    } else {
      for (int16 tip = 0; tip != n_tip; tip++) {
        t2_descendantsof[index2 * n_tip + tip] |= 
          t2_descendantsof[(child2 - n_tip) * n_tip + tip];
      }
    }
  }
  
  int_fast16_t *M = new int_fast16_t[n_all_nodes * n_all_nodes]();
  for (int16 i = 0; i != n_edge; i++) {
    const int16 node1 = edge1(i, 1);
    for (int16 j = 0; j != n_edge; j++) {
      fill_M(node1, edge2(j, 1), M, t1_descendantsof,
             t2_descendantsof, t1_left, t1_right,
             t2_left, t2_right, n_tip, n_all_nodes);
    }
  }
  // Root node:
  fill_M(n_tip, n_tip, M, t1_descendantsof,
         t2_descendantsof, t1_left, t1_right,
         t2_left, t2_right, n_tip, n_all_nodes);
  
  delete[] t1_descendantsof;
  delete[] t2_descendantsof;
  int ret = int(M[n_tip * n_all_nodes + n_tip]);
  delete[] M;
  
  return ret;
}
