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
             int16 M[][MAST_MAX_ALLNODE],
             const bool t1_descendantsof[][MAST_MAX_TIP],
             const bool t2_descendantsof[][MAST_MAX_TIP],
             const int16* t1_left,
             const int16* t1_right,
             const int16* t2_left,
             const int16* t2_right,
             const int16 n_tip) {
  M[node1][node2] = 0; // TODO remove this line and rewrite all subsequent assignemnts as ?1:0
  if (node1 < n_tip) { // Node 1 is a leaf
    if (node2 < n_tip) { // Node 1 and node 2 are leaves
      if (node1 == node2) {
        M[node1][node2] = 1;
      }
    } else { // Node 1 is a leaf; node 2 is internal
      if (t2_descendantsof[node2 - n_tip][node1]) {
        M[node1][node2] = 1;
      }
    }
  } else {
    if (node2 < n_tip) { // Node 1 is internal; node 2 is a leaf
      if (t1_descendantsof[node1 - n_tip][node2]) {
        M[node1][node2] = 1;
      }
    } else { // Node 1 & Node 2 are internal
      const int16 l1 = t1_left[node1 - n_tip],
                  r1 = t1_right[node1 - n_tip],
                  l2 = t2_left[node2 - n_tip],
                  r2 = t2_right[node2 - n_tip];
      M[node1][node2] = max_of_six(
        M[l1][l2] + M[r1][r2],
        M[l1][r2] + M[r1][l2],
        M[node1][l2],
        M[node1][r2],
        M[l1][node2],
        M[r1][node2]);
    }
  }
}


// edge1 and edge2 should be edge matrices taken from trees of class `phylo`,
// with one subtracted from all entries so that leaves are numbered from 
// 0..(nTip - 1).
// [[Rcpp::export]]
int cpp_mast (IntegerMatrix edge1, IntegerMatrix edge2, IntegerVector nTip) {
  const int16 n_tip = nTip[0], n_edge = edge1.nrow();
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
  bool t1_descendantsof[MAST_MAX_NODE][MAST_MAX_TIP] = {},
       t2_descendantsof[MAST_MAX_NODE][MAST_MAX_TIP] = {};
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
      t1_descendantsof[index1][child1] = true;
    } else { // Child 1 is a node
      for (int16 tip = 0; tip != n_tip; tip++) {
        t1_descendantsof[index1][tip] |= 
          t1_descendantsof[child1 - n_tip][tip];
      }
    }
    if (child2 < n_tip) {
      t2_descendantsof[index2][child2] = true;
    } else {
      for (int16 tip = 0; tip != n_tip; tip++) {
        t2_descendantsof[index2][tip] |= 
          t2_descendantsof[child2 - n_tip][tip];
      }
    }
  }
  
  int16 M[MAST_MAX_ALLNODE][MAST_MAX_ALLNODE];
  for (int16 i = 0; i != n_edge; i++) {
    int16 node1 = edge1(i, 1);
    for (int16 j = 0; j != n_edge; j++) {
      fill_M(node1, edge2(j, 1), M, t1_descendantsof,
             t2_descendantsof, t1_left, t1_right,
             t2_left, t2_right, n_tip);
    }
  }
  // Root node:
  fill_M(n_tip, n_tip, M, t1_descendantsof,
         t2_descendantsof, t1_left, t1_right,
         t2_left, t2_right, n_tip);
  
  return int(M[n_tip][n_tip]);
}
