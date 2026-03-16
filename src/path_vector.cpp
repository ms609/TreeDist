#include <Rcpp/Lightest>
#include <TreeTools/renumber_tree.h> // for postorder_order
#include <algorithm> // for std::copy
#include <memory> // for make_unique
using namespace Rcpp;

#define PO_PARENT(i) edge(postorder[i] - 1, 0)
#define PO_CHILD(i) edge(postorder[i] - 1, 1)

#define ANC(label, i) ancestry[n_tip * (label - 1) + i]

// [[Rcpp::export]]
IntegerVector path_vector(IntegerMatrix edge) {
  
  IntegerVector postorder = TreeTools::postorder_order(edge);
  
  const int n_edge = edge.nrow();
  const int n_vert = n_edge + 1;
  const int root_node = PO_PARENT(n_edge - 1);
  const int n_tip = root_node - 1;

  auto ancestry = std::make_unique<int[]>(n_tip * n_vert);
  auto n_ancs = std::make_unique<int[]>(n_vert + 1);
  
  for (int i = n_edge; i--; ) { // Preorder traversal
    const int parent = PO_PARENT(i);
    const int child = PO_CHILD(i);
    const int parent_ancs = n_ancs[parent];
    
    const int* anc_parent = &ancestry[n_tip * (parent - 1)];
    int* anc_child = &ancestry[n_tip * (child - 1)];
    
    anc_child[parent_ancs] = child;
    n_ancs[child] = parent_ancs + 1;
    
    std::copy(anc_parent, anc_parent + parent_ancs, anc_child);
  }
  
  const int ret_size = n_tip * (n_tip - 1) / 2;
  int ptr = ret_size;
  IntegerVector ret(ptr);
  for (int i = n_tip - 1; i--; ) { // faster in reverse, for a benchmark
    for (int j = n_tip - i - 1; j--; ) {
      const int tip_i = i + 1;
      const int tip_j = i + 2 + j;
      const int ancs_i = n_ancs[tip_i];
      const int ancs_j = n_ancs[tip_j];
      const int min_ancs = ancs_i > ancs_j ? ancs_j : ancs_i;
      
      assert(ptr > 0);
      assert(tip_j > tip_i);
      
      const int* anc_i_base = &ancestry[n_tip * (tip_i - 1)];
      const int* anc_j_base = &ancestry[n_tip * (tip_j - 1)];
      
      int common = 0;
      for (; common != min_ancs; ++common) {
        if (anc_i_base[common] != anc_j_base[common]) {
          break;
        }
      }
      --ptr;
      ret[ptr] = n_ancs[tip_i] + n_ancs[tip_j] - (common << 1);
    }
  }
  
  return ret;
}
