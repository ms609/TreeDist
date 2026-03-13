#include <Rcpp/Lightest>
#include <TreeTools/renumber_tree.h> // for postorder_order
#include <algorithm> // for std::copy
#include <cmath> // for sqrt
#include <memory> // for make_unique
#ifdef _OPENMP
#include <omp.h>
#endif
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

// Kendall-Colijn vector: for each pair of tips (in tip_order), return the
// depth of their MRCA from the root (= number of shared prefix entries in
// the ancestry arrays).
// [[Rcpp::export]]
IntegerVector cpp_kc_vector(IntegerMatrix edge, IntegerVector tip_order) {

  IntegerVector postorder = TreeTools::postorder_order(edge);

  const int n_edge = edge.nrow();
  const int n_vert = n_edge + 1;
  const int root_node = PO_PARENT(n_edge - 1);
  const int n_tip = root_node - 1;

  // Build ancestry matrix (same as path_vector)
  auto ancestry = std::make_unique<int[]>(n_tip * n_vert);
  auto n_ancs = std::make_unique<int[]>(n_vert + 1);

  for (int i = n_edge; i--; ) {
    const int parent = PO_PARENT(i);
    const int child = PO_CHILD(i);
    const int parent_ancs = n_ancs[parent];

    const int* anc_parent = &ancestry[n_tip * (parent - 1)];
    int* anc_child = &ancestry[n_tip * (child - 1)];

    anc_child[parent_ancs] = child;
    n_ancs[child] = parent_ancs + 1;

    std::copy(anc_parent, anc_parent + parent_ancs, anc_child);
  }

  // Iterate pairs in tip_order (matching R's combn(tipOrder, 2) column order)
  const int ret_size = n_tip * (n_tip - 1) / 2;
  IntegerVector ret(ret_size);
  int ptr = 0;
  for (int i = 0; i < n_tip - 1; ++i) {
    const int tip_i = tip_order[i];
    const int ancs_i = n_ancs[tip_i];
    const int* anc_i_base = &ancestry[n_tip * (tip_i - 1)];

    for (int j = i + 1; j < n_tip; ++j) {
      const int tip_j = tip_order[j];
      const int ancs_j = n_ancs[tip_j];
      const int min_ancs = ancs_i > ancs_j ? ancs_j : ancs_i;

      const int* anc_j_base = &ancestry[n_tip * (tip_j - 1)];

      int common = 0;
      for (; common != min_ancs; ++common) {
        if (anc_i_base[common] != anc_j_base[common]) {
          break;
        }
      }
      ret[ptr++] = common;
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
  const int n_pairs = col1 * col2;
  const int* v1 = INTEGER(vec1);
  const int* v2 = INTEGER(vec2);
  double* out = REAL(ret);

  #pragma omp parallel for schedule(dynamic) if(n_pairs > 100)
  for (int idx = 0; idx < n_pairs; ++idx) {
    const int i = idx / col2;
    const int j = idx % col2;
    const int* col_i = v1 + (std::size_t)i * n_row;
    const int* col_j = v2 + (std::size_t)j * n_row;
    int64_t val = 0;
    for (int row = 0; row < n_row; ++row) {
      const int64_t x = col_i[row] - col_j[row];
      val += x * x;
    }
    out[(std::size_t)j * col1 + i] = std::sqrt((double)val);
  }

  return ret;
}

// dist-order index for pair (i, j) where i < j, in an n×n lower-triangle:
//   offset = n*i - i*(i+1)/2 + (j - i - 1)
static inline int dist_index(int i, int j, int n) {
  return n * i - i * (i + 1) / 2 + (j - i - 1);
}

// [[Rcpp::export]]
NumericVector pair_diff_euclidean(const IntegerMatrix vecs) {
  const int
    n_col = vecs.cols(),
    n_row = vecs.rows()
  ;

  const int n_pairs = n_col * (n_col - 1) / 2;
  NumericVector ret(n_pairs);
  const int* v = INTEGER(vecs);
  double* out = REAL(ret);

  // Linearise the upper triangle: pair index k maps to (i, j) with i < j.
  // Parallelise over the outer index i (each i owns a contiguous block of
  // output slots), avoiding the need to map linear k → (i,j).
  #pragma omp parallel for schedule(dynamic) if(n_pairs > 100)
  for (int i = 0; i < n_col - 1; ++i) {
    const int* col_i = v + (std::size_t)i * n_row;
    for (int j = i + 1; j < n_col; ++j) {
      const int* col_j = v + (std::size_t)j * n_row;
      int64_t val = 0;
      for (int row = 0; row < n_row; ++row) {
        const int64_t x = col_i[row] - col_j[row];
        val += x * x;
      }
      out[dist_index(i, j, n_col)] = std::sqrt((double)val);
    }
  }

  return ret;
}
