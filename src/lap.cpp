/************************************************************************
 *
 *  lap.cpp — R entry point (lapjv) for the Linear Assignment Problem.
 *
 *  The lap() algorithm itself lives in inst/include/TreeDist/lap_impl.h,
 *  the single source of truth shared with downstream LinkingTo packages.
 *  See that header for the Jonker–Volgenant attribution and provenance.
 *
 *************************************************************************/

#include "lap.h"
#include <Rcpp/Lightest>

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
  
  // Build cost matrix.  Fill the transposed buffer first (matching R's
  // column-major storage for sequential reads) then untranspose.
  cost_matrix input(max_dim);
  const double* __restrict__ src_data = REAL(x);
  cost* __restrict__ t_ptr = input.col(0);
  const std::size_t dim8 = input.dim8();
  
  for (lap_col c = 0; c < n_col; ++c) {
    const std::size_t t_off = static_cast<std::size_t>(c) * dim8;
    const std::size_t s_off = static_cast<std::size_t>(c) * n_row;
    for (lap_row r = 0; r < n_row; ++r) {
      t_ptr[t_off + r] = static_cast<cost>(src_data[s_off + r] * scale_factor);
    }
    // Pad remaining rows in this transposed column
    for (lap_row r = n_row; r < max_dim; ++r) {
      t_ptr[t_off + r] = max_score;
    }
    for (std::size_t r = max_dim; r < dim8; ++r) {
      t_ptr[t_off + r] = max_score;
    }
  }
  // Pad remaining transposed columns
  for (lap_col c = n_col; c < max_dim; ++c) {
    const std::size_t t_off = static_cast<std::size_t>(c) * dim8;
    for (std::size_t r = 0; r < dim8; ++r) {
      t_ptr[t_off + r] = max_score;
    }
  }
  
  // Untranspose: t_data_ -> data_
  input.makeUntranspose();
  
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

// Pull in the one true definition of TreeDist::lap().  Wire the interrupt
// macro to Rcpp so the package keeps its user-interrupt checks.
#define TREEDIST_CHECK_INTERRUPT() Rcpp::checkUserInterrupt()
#define TREEDIST_LAP_IMPLEMENTATION
#include <TreeDist/lap_impl.h>
