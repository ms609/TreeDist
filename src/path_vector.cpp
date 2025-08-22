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
  const int
    n_edge = edge.nrow(),
    n_vert = n_edge + 1
  ;
  IntegerVector postorder = TreeTools::postorder_order(edge);
  const int
    root_node = PO_PARENT(n_edge - 1),
    n_tip = root_node - 1
  ;

  auto ancestry = std::make_unique<int[]>(n_tip * n_vert);
  auto n_ancs = std::make_unique<int[]>(n_vert + 1);
  
  for (int i = n_edge; i--; ) { // Preorder traversal
    const int
      parent = PO_PARENT(i),
      child = PO_CHILD(i)
    ;
    ANC(child, n_ancs[parent]) = child;
    n_ancs[child] = n_ancs[parent] + 1;
    for (int j = n_ancs[parent]; j--; ) {
      ANC(child, j) = ANC(parent, j);
    }
  }
  
  int ptr = n_tip * (n_tip - 1) / 2;
  IntegerVector ret(ptr);
  for (int i = n_tip - 1; i--; ) {
    for (int j = n_tip - i - 1; j--; ) {
      const int 
        tip_i = i + 1,
        tip_j = i + 2 + j,
        ancs_i = n_ancs[tip_i],
        ancs_j = n_ancs[tip_j],
        min_ancs = ancs_i > ancs_j ? ancs_j : ancs_i
      ;
      
      assert(ptr > 0);
      assert(tip_j > tip_i);
      
      // Rcout << tip_i << " (" << n_ancs[tip_i] << "), "
      //       << tip_j << " (" << n_ancs[tip_j] << "): ";
      int common = 0;
      for (; common != min_ancs; ++common) {
        // Rcout << "{" << ANC(tip_i, common) <<" " << ANC(tip_j, common) << "} ";
        if (ANC(tip_i, common) != ANC(tip_j, common)) {
          // Rcout << ANC(tip_i, common) << " != " << ANC(tip_j, common)
          //       << "; common = " << common << ";\n";
          break;
        }
      }
      ret[--ptr] = n_ancs[tip_i] + n_ancs[tip_j] - common - common;
    }
  }
  
  return ret;
}

// [[Rcpp::export]]
NumericMatrix vec_diff_euclidean(const IntegerMatrix vec1,
                                 const IntegerMatrix vec2) {
  const int 
    col1 = vec1.cols(),
    col2 = vec2.cols(),
    n_row = vec1.rows()
  ;
  assert(n_row == vec2.rows());
  
  NumericMatrix ret(col1, col2);
  for (int i = col1; i--; ) {
    for (int j = col2; j--; ) {
      int val = 0;
      for (int row = n_row; row--; ) {
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
