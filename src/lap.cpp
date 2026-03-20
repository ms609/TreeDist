/************************************************************************
 *
 *  lap.cpp
 version 1.0 - 4 September 1996
 author: Roy Jonker @ MagicLogic Optimization Inc.
 e-mail: roy_jonker@magiclogic.com
 
 Code for Linear Assignment Problem, according to
 
 "A Shortest Augmenting Path Algorithm for Dense and Sparse Linear
 Assignment Problems," Computing 38, 325-340, 1987
 
 by
 
 R. Jonker and A. Volgenant, University of Amsterdam.
 
 *
 CHANGED 2025-08-01 by Martin Smith <martin.smith@durham.ac.uk> to optimize
 performance.
 CHANGED 2020-01-01 by Martin Smith <martin.smith@durham.ac.uk> for 
 integration with R.
 CHANGED 2016-05-13 by Yong Yang(yongyanglink@gmail.com) in column reduction 
 part according to matlab version of LAPJV algorithm
 https://github.com/yongyanghz/LAPJV-algorithm-c/blob/master/src/lap.cpp
 (Copyright (c) 2010, Yi Cao All rights reserved)--
 https://www.mathworks.com/matlabcentral/fileexchange/26836-lapjv-jonker-volgenant-algorithm-for-linear-assignment-problem-v3-0:
 *
 *************************************************************************/

#include <Rcpp/Lightest>

// Provide the LAP implementation in this translation unit.
#define TREEDIST_CHECK_INTERRUPT() Rcpp::checkUserInterrupt()
#define TREEDIST_LAP_IMPLEMENTATION
#include <TreeDist/lap_impl.h>

#include "lap.h"

// [[Rcpp::export]]
Rcpp::List lapjv(Rcpp::NumericMatrix &x, Rcpp::NumericVector &maxX) {
  const lap_dim n_row = x.nrow();
  const lap_dim n_col = x.ncol();
  const lap_dim max_dim = std::max(n_row, n_col);
  const lap_dim spare_rows = n_row - n_col;
  const cost max_score = cost(BIG / max_dim);
  const double x_max = maxX[0];
  const double scale_factor = max_score / x_max;
  
  std::vector<lap_col> rowsol(max_dim);
  std::vector<lap_row> colsol(max_dim);
  
  // Build cost matrix from R's column-major NumericMatrix.
  cost_matrix input(max_dim);
  const double* __restrict__ src_data = REAL(x);
  for (lap_row r = 0; r < n_row; ++r) {
    for (lap_col c = 0; c < n_col; ++c) {
      input(r, c) = static_cast<cost>(
        src_data[static_cast<std::size_t>(c) * n_row + r] * scale_factor);
    }
    input.padRowAfterCol(r, n_col, max_score);
  }
  input.padAfterRow(n_row, max_score);
  
  cost score = lap(max_dim, input, rowsol, colsol);
  
  std::vector<int> matching;
  matching.reserve(n_row);
  
  for (lap_dim i = 0; i < n_row; ++i) {
    const int match = (rowsol[i] < n_col) ? rowsol[i] + 1 : NA_INTEGER;
    matching.push_back(match);
  }
  
  return Rcpp::List::create(
    Rcpp::Named("score") = (static_cast<double>(score) -
      (std::abs(spare_rows) * max_score)) / max_score * x_max,
    Rcpp::_["matching"] = matching
  );
}
