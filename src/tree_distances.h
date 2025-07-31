#ifndef _TREEDIST_TREE_DISTANCES_H
#define _TREEDIST_TREE_DISTANCES_H

#include <TreeTools/SplitList.h>
#include <Rcpp/Lightest>

#include <limits> /* for numeric_limits */
#include <vector> /* for vector */

#include "ints.h"

/*************** TYPES      *******************/

using cost = int_fast64_t;
using lap_dim = int16;
using lap_row = lap_dim;
using lap_col = lap_dim;

/***** Constants requiring initialization *****/

constexpr splitbit ALL_ONES = (std::numeric_limits<splitbit>::max)();
extern double lg2[int32(SL_MAX_TIPS - 1) * (SL_MAX_TIPS - 1) + 1];
extern double lg2_double_factorial[SL_MAX_TIPS + SL_MAX_TIPS - 2];
extern double lg2_unrooted[SL_MAX_TIPS + 2];
extern double *lg2_rooted;

/* For a reason I've not determined, shrinking BIG is necessary to avoid 
 * an infinite loop in lap. */
constexpr cost BIG = (std::numeric_limits<cost>::max)() / SL_MAX_SPLITS;
constexpr cost ROUND_PRECISION = 2048 * 2048;

template<typename T>
class FlatMatrix {
private:
  static constexpr size_t BLOCK_SIZE = 8;
  const size_t dim_; // Important not to use int16, which will overflow on *
  const size_t dim8_;
  alignas(64) std::vector<T> data_;
  const size_t block_containing(const size_t x) {
    return ((x + BLOCK_SIZE - 1) / BLOCK_SIZE) * BLOCK_SIZE;
  }
  
public:
  FlatMatrix(size_t dim)
    : dim_(dim),
      dim8_(block_containing(dim_)),
      data_(std::vector<T>(dim8_ * dim8_)) {}
      
  FlatMatrix(const Rcpp::NumericMatrix& src, const double x_max)
    : dim_((std::max(src.nrow(), src.ncol()))),  // or pad here as needed
      dim8_(block_containing(dim_)),
      data_(std::vector<T>(dim8_ * dim_))
  {
    // Compute scale factor
    const cost max_score = cost(BIG / dim_);
    double scale_factor = max_score / x_max;
    
    const int16 n_row = src.nrow();
    const int16 n_col = src.ncol();
    
    for (int16 r = 0; r < n_row; ++r) {
      const size_t r_offset = static_cast<size_t>(r) * dim8_;
      
      for (int16 c = 0; c < n_col; ++c) {
        data_[r_offset + c] = static_cast<T>(src(r, c) * scale_factor);
      }
      
      padRowAfterCol(r, n_col, max_score);
    }
    
    padAfterRow(n_row, max_score);
  }
  
  // Access operator for read/write
  T& operator()(lap_row i, lap_col j) {
    return data_[static_cast<size_t>(i) * dim8_ + j];
  }
  
  // Const version for read-only access
  [[nodiscard]] const T& operator()(lap_row i, lap_col j) const {
    return data_[static_cast<size_t>(i) * dim8_ + j];
  }
  
  [[nodiscard]] const T& row0(lap_col j) const {
    return data_[j];
  }
  
  [[nodiscard]] const T& entry0(lap_row i) const {
    return data_[static_cast<size_t>(i) * dim8_];
  }
  
  [[nodiscard]] T* row(lap_row i) {
    return &data_[static_cast<size_t>(i) * dim8_];
  }
  
  [[nodiscard]] const T* row(lap_row i) const {
    return &data_[static_cast<size_t>(i) * dim8_];
  }
  
  void padAfterRow(lap_row start_row, T value) {
    size_t start_index = static_cast<size_t>(start_row) * dim8_;
    std::fill(data_.begin() + start_index, data_.end(), value);
  }
  
  void padRowAfterCol(const lap_row r, const lap_col start_col,
                      const T value) {
    size_t r_offset = r * dim8_;
    size_t actual_start_col = static_cast<size_t>(start_col);
    size_t start_index = r_offset + actual_start_col;
    size_t end_index = start_index + dim_ - actual_start_col;
    std::fill(data_.begin() + start_index, data_.begin() + end_index, value);
  }
};

using cost_matrix = FlatMatrix<cost>;

/*************** FUNCTIONS  *******************/

extern cost lap(int16 dim,
                cost_matrix &input_cost,
                std::vector<lap_col> &rowsol,
                std::vector<lap_row> &colsol);

extern double 
  mmsi_score(const int16 n_same, const int16 n_a_and_b,
             const int16 n_different, const int16 n_a_only),
  ic_element(const int16 nkK, const int16 nk,
             const int16 nK, const int16 n),
  ic_matching(const int16 a, const int16 b, const int16 n),
  one_overlap(const int16 a, const int16 b, const int16 n),
  one_overlap_notb(const int16 a, const int16 n_minus_b, const int16 n),
  spi_overlap(const splitbit* a_state, const splitbit* b_state,
              const int16 n_tips, const int16 in_a, const int16 in_b,
              const int16 n_bins);

extern Rcpp::List cpp_robinson_foulds_distance(Rcpp::RawMatrix x,
                                               Rcpp::RawMatrix y,
                                               Rcpp::IntegerVector nTip);
#endif