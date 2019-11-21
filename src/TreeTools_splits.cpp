#include "splits.h"

const int powers_of_two[16] = {1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024,
                               2048, 4096, 8192, 16384, 32768};
const int BIN_SIZE = 8;

// Edges must be listed in postorder
// [[Rcpp::export]]
RawMatrix cpp_edge_to_splits(IntegerMatrix edge, IntegerVector nTip) {
  if (edge.cols() != 2) {
    throw std::invalid_argument("Edge matrix must contain two columns");
  }
  
  const int n_edge = edge.rows(),
    n_node = n_edge + 1,
    n_tip = nTip[0],
                n_return = n_edge - n_tip - 1,
                n_bin = ((n_tip - 1) / BIN_SIZE) + 1;
  
  if (n_edge == n_tip) { /* No internal nodes resolved */
    return RawMatrix (0, n_bin);
  }
  
  int** splits = new int*[n_node];
  for (int i = 0; i != n_node; i++) {
    splits[i] = new int[n_bin];
    for (int j = 0; j != n_bin; j++) {
      splits[i][j] = 0;
    }
  }
  
  for (int i = 0; i != n_tip; i++) {
    splits[i][(int) i / BIN_SIZE] = powers_of_two[i % BIN_SIZE];
  }
  
  for (int i = 0; i != n_edge - 1; i++) { /* final edge is second root edge */
    for (int j = 0; j != n_bin; j++) {
      splits[(int) edge(i, 0) - 1][j] |= splits[(int) edge(i, 1) - 1][j];
    }
  }
  
  RawMatrix ret(n_return, n_bin);
  IntegerVector names(n_return);
  int n_trivial = 0;
  const int trivial_one = edge(n_edge - 1, 0) - 1,
    trivial_two = ((edge(n_edge - 1, 1) <= n_tip) ?
                     edge(n_edge - 2, 1) : edge(n_edge - 1, 1)) - 1;
  for (int i = n_tip; i != n_node; i++) {
    if (i == trivial_one || i == trivial_two) {
      n_trivial++;
    } else {
      for (int j = 0; j != n_bin; j++) {
        ret(i - n_tip - n_trivial, j) = splits[i][j];
        names[i - n_tip - n_trivial] = (i + 1);
      }
    }
  }
  
  rownames(ret) = names;
  return(ret);
}
