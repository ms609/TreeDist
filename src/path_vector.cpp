#include <Rcpp/Lightest>
#include <TreeTools/renumber_tree.h> // for postorder_order
#include <cmath> // for sqrt
#include <memory> // for make_unique
using namespace Rcpp;

#define PO_PARENT(i) edge(postorder[i] - 1, 0)
#define PO_CHILD(i) edge(postorder[i] - 1, 1)
#define GET_DIST(i, j) i < j ? dist_from(j, i) : dist_from(i, j)
#define SET_DIST(i, j, x) i < j ? dist_from(j, i) = x : dist_from(i, j) = x

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
    
    for (int j = 0; j < parent_ancs; ++j) {
      anc_child[j] = anc_parent[j];
    }
  }
  
  int ptr = n_tip * (n_tip - 1) / 2;
  IntegerVector ret(ptr);
  for (int i = n_tip - 1; i--; ) {
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

// [[Rcpp::export]]
NumericMatrix vec_diff_euclidean(const IntegerMatrix vec1,
                                 const IntegerMatrix vec2) {
  const int col1 = vec1.cols();
  const int col2 = vec2.cols();
  const int n_row = vec1.rows();
  
  assert(n_row == vec2.rows());
  
  NumericMatrix ret(col1, col2);
  for (int i = 0; i < col1; ++i) {
    for (int j = 0; j < col2; ++j) {
      int val = 0;
      for (int row = 0; row < n_row; ++row) {
        const int x = vec1(row, i) - vec2(row, j);
        val += x * x;
      }
      ret(i, j) = std::sqrt(val);
    }
  }
  
  return ret;
}

// [[Rcpp::export]]
NumericVector pair_diff_euclidean(const IntegerMatrix vecs) {
  const int 
    n_col = vecs.cols(),
    n_row = vecs.rows()
  ;
  
  int ptr = n_col * (n_col - 1) / 2;
  NumericVector ret(ptr);
  for (int i = n_col - 1; i--; ) {
    for (int j_it = n_col - i - 1; j_it--; ) {
      const int j = i + 1 + j_it;
      assert(ptr > 0);
      assert(j > i);
      
      int val = 0;
      for (int row = n_row; row--; ) {
        const int x = vecs(row, i) - vecs(row, j);
        val += x * x;
      }
      ret[--ptr] = std::sqrt(val);
    }
  }
  
  return ret;
}
