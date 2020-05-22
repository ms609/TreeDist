#include <Rcpp.h>
#include "ints.h"
using namespace Rcpp;

const uint16 MAX_MAST_NODE = 1023, /* 4095 causes crash */
  MAX_MAST_TIP = MAX_MAST_NODE + 1,
  MAX_MAST_ALLNODE = (MAX_MAST_NODE * 2) - 1;

uint16 max_of_six (const uint16 m1,
                   const uint16 m2,
                   const uint16 m3,
                   const uint16 m4,
                   const uint16 m5,
                   const uint16 m6) {
  uint16 largest = (m1 > m2 ? m1 : m2);
  if (m3 > largest) largest = m3;
  if (m4 > largest) largest = m4;
  if (m5 > largest) largest = m5;
  return (m6 > largest ? m6 : largest);
}

void fill_M (const uint16 node1, const uint16 node2,
             uint16 M[][MAX_MAST_NODE],
             const bool t1_descendantsof[][MAX_MAST_TIP],
             const bool t2_descendantsof[][MAX_MAST_TIP],
             const uint16* t1_left,
             const uint16* t1_right,
             const uint16* t2_left,
             const uint16* t2_right,
             const uint16 n_tip) {
  M[node1][node2] = 0;
  if (node1 < n_tip) {
    if (node2 < n_tip) {
      if (node1 == node2) {
        M[node1][node2] = 1;
      }
    } else {
      if (t2_descendantsof[node2 - n_tip][node1]) {
        M[node1][node2] = 1;
      }
    }
  } else {
    if (node2 < n_tip) {
      if (t1_descendantsof[node1 - n_tip][node2]) {
        M[node1][node2] = 1;
      }
    } else {
      const uint16 l1 = t1_left[node1],
                              r1 = t1_right[node1],
                              l2 = t2_left[node2],
                              r2 = t2_right[node2];
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
  const uint16 n_tip = nTip[0], n_edge = edge1.nrow();
  uint16 
    t1_left[MAX_MAST_NODE] = {}, t1_right[MAX_MAST_NODE] = {},
    t2_left[MAX_MAST_NODE] = {}, t2_right[MAX_MAST_NODE] = {};
  bool t1_has_child[MAX_MAST_NODE] = {}, t2_has_child[MAX_MAST_NODE] = {};
  bool t1_descendantsof[MAX_MAST_NODE][MAX_MAST_TIP] = {},
       t2_descendantsof[MAX_MAST_NODE][MAX_MAST_TIP] = {};
  if (edge2.nrow() != int16(n_edge)) {
    throw std::length_error("Both trees must contain the same number of edges.");
  }
  for (uint16 i = 0; i < n_edge; i++) {
    const uint16 
      parent1 = edge1(i, 0), child1 = edge1(i, 1),
      parent2 = edge2(i, 0), child2 = edge2(i, 1);
    if (t1_has_child[parent1]) {
      t1_right[parent1] = child1;
    } else {
      t1_left[parent1] = child1;
      t1_has_child[parent1] = true;
    }
    if (t2_has_child[parent2]) {
      t2_right[parent2] = child2; 
    } else {
      t2_left[parent2] = child2;
      t2_has_child[parent2] = true;
    }
    
    if (child1 < n_tip) {
      t1_descendantsof[parent1 - n_tip][child1] = true;
    } else {
      for (uint16 tip = 0; tip < n_tip; tip++) {
        t1_descendantsof[parent1 - n_tip][tip] |= 
          t1_descendantsof[child1 - n_tip][tip];
      }
    }
    if (child2 < n_tip) {
      t2_descendantsof[parent2 - n_tip][child2] = true;
    } else {
      for (uint16 tip = 0; tip < n_tip; tip++) {
        t2_descendantsof[parent2 - n_tip][tip] |= 
          t2_descendantsof[child2 - n_tip][tip];
      }
    }
  }
  
  uint16 M[MAX_MAST_NODE][MAX_MAST_NODE];
  for (uint16 i = 0; i < n_edge; i++) {
    uint16 node1 = edge1(i, 1);
    for (uint16 j = 0; j < n_edge; j++) {
      fill_M(node1, edge2(j, 1), M, t1_descendantsof,
             t2_descendantsof, t1_left, t1_right,
             t2_left, t2_right, n_tip);
    }
  }
  fill_M(n_tip, n_tip, M, t1_descendantsof,
         t2_descendantsof, t1_left, t1_right,
         t2_left, t2_right, n_tip);
  
  return int(M[n_tip][n_tip]);
}
