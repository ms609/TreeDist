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
using row_offset = size_t; // must hold int16 * int16

template<typename T>
class FlatMatrix {
private:
  size_t dim_; // Important not to use int16, which will overflow
  alignas(64) std::vector<T> data_;
  
public:
  FlatMatrix(size_t dim)
    : data_(std::vector<T>(dim * dim)), dim_(dim) {}
  
  // Access operator for read/write
  T& operator()(lap_row i, lap_col j) {
    return data_[static_cast<size_t>(i) * dim_ + j];
  }
  
  // Const version for read-only access
  [[nodiscard]] const T& operator()(lap_row i, lap_col j) const {
    return data_[static_cast<size_t>(i) * dim_ + j];
  }
  
  [[nodiscard]] const T& row0(lap_col j) const {
    return data_[j];
  }
  
  // Access operator for read/write
  [[nodiscard]] T& fromRow(row_offset i, lap_col j) {
    return data_[i + j];
  }
  
  // Const version for read-only access
  [[nodiscard]] const T& fromRow(row_offset i, lap_col j) const {
    return data_[i + j];
  }
  
  [[nodiscard]] const T& entry0(lap_row i) const {
    return data_[static_cast<size_t>(i) * dim_];
  }
  
  [[nodiscard]] const T& atOffset(row_offset i) const {
    return data_[i];
  }
  
  [[nodiscard]] T* row(lap_row i) {
    return &data_[static_cast<size_t>(i) * dim_];
  }
  
  [[nodiscard]] const T* row(lap_row i) const {
    return &data_[static_cast<size_t>(i) * dim_];
  }
  
  void padAfterRow(lap_row start_row, T value) {
    size_t start_index = static_cast<size_t>(start_row) * dim_;
    std::fill(data_.begin() + start_index, data_.end(), value);
  }
  
  void padRowAfterCol(const row_offset r_offset, const lap_col start_col,
                      const T value) {
    size_t actual_start_col = static_cast<size_t>(start_col);
    size_t start_index = r_offset + actual_start_col;
    size_t end_index = start_index + dim_ - actual_start_col;
    std::fill(data_.begin() + start_index, data_.begin() + end_index, value);
  }
};

using cost_matrix = FlatMatrix<cost>;

constexpr splitbit ALL_ONES = (std::numeric_limits<splitbit>::max)();

/* For a reason I've not determined, shrinking BIG is necessary to avoid 
 * an infinite loop in lap. */
constexpr cost BIG = (std::numeric_limits<cost>::max)() / SL_MAX_SPLITS;
constexpr cost ROUND_PRECISION = 2048 * 2048;

/***** Constants requiring initialization *****/

extern double lg2[int32(SL_MAX_TIPS - 1) * (SL_MAX_TIPS - 1) + 1];
extern double lg2_double_factorial[SL_MAX_TIPS + SL_MAX_TIPS - 2];
extern double lg2_unrooted[SL_MAX_TIPS + 2];
extern double *lg2_rooted;

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