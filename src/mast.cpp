#include <Rcpp.h>
using namespace Rcpp;

const unsigned int MAX_MAST_NODE = 4095,
  MAX_MAST_TIP = MAX_MAST_NODE + 1,
  MAX_MAST_ALLNODE = MAX_MAST_NODE * 2 - 1;

unsigned int max_of_six (const unsigned int m1,
                         const unsigned int m2,
                         const unsigned int m3,
                         const unsigned int m4,
                         const unsigned int m5,
                         const unsigned int m6) {
  unsigned int largest = (m1 > m2 ? m1 : m2);
  if (m3 > largest) largest = m3;
  if (m4 > largest) largest = m4;
  if (m5 > largest) largest = m5;
  return (m6 > largest ? m6 : largest);
}

// [[Rcpp::export]]
int cpp_mast (IntegerMatrix edge1, IntegerMatrix edge2, IntegerVector nTip) {
  const unsigned int n_tip = nTip[0], n_edge = edge1.size(), n_node = n_tip - 1,
    all_node = n_tip + n_node;
  unsigned int 
    t1_left[MAX_MAST_NODE] = {}, t1_right[MAX_MAST_NODE] = {},
    t2_left[MAX_MAST_NODE] = {}, t2_right[MAX_MAST_NODE] = {};
  bool t1_has_child[MAX_MAST_NODE] = {}, t2_has_child[MAX_MAST_NODE] = {};
  bool t1_descendantsof[MAX_MAST_NODE][MAX_MAST_TIP] = {},
       t2_descendantsof[MAX_MAST_NODE][MAX_MAST_TIP] = {};
  
  if (edge2.size() != n_edge) {
    throw std::length_error("Both trees must contain the same number of edges.");
  }
  
  for (unsigned int i = 0; i < n_edge; i++) {
    const unsigned int parent1 = edge1(i, 0), child1 = edge1(i, 1),
      parent2 = edge2(i, 0), child2 = edge2(i, 1);
    if (t1_has_child[parent1]) {
          if (t1_descendantsof[node1][node2]) {
            M[node1][node2] = 1;
          }
        } else {
          const unsigned int l1 = t1_left[node1],
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
  }
  return (M[n_tip + 1][n_tip + 1]);
}
