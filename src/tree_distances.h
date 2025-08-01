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
constexpr size_t BLOCK_SIZE = 16;

class CostMatrix {
private:
  const size_t dim_; // Important not to use int16, which will overflow on *
  const size_t dim8_;
  alignas(64) std::vector<cost> data_;
  alignas(64) std::vector<cost> t_data_;
  const size_t block_containing(const size_t x) {
    return ((x + BLOCK_SIZE - 1) / BLOCK_SIZE) * BLOCK_SIZE;
  }
  
public:
  CostMatrix(size_t dim)
    : dim_(dim),
      dim8_(block_containing(dim_)),
      data_(std::vector<cost>(dim8_ * dim8_)),
      t_data_(std::vector<cost>(dim8_ * dim8_)) {}
      
  CostMatrix(const Rcpp::NumericMatrix& src, const double x_max)
    : dim_((std::max(src.nrow(), src.ncol()))),  // or pad here as needed
      dim8_(block_containing(dim_)),
      data_(std::vector<cost>(dim8_ * dim_)),
      t_data_(std::vector<cost>(dim8_ * dim_))
  {
    // Compute scale factor
    const cost max_score = cost(BIG / dim_);
    double scale_factor = max_score / x_max;
    
    const lap_row n_row = src.nrow();
    const lap_col n_col = src.ncol();
    
    for (lap_row r = 0; r < n_row; ++r) {
      const size_t r_offset = static_cast<size_t>(r) * dim8_;
      
      for (lap_col c = 0; c < n_col; ++c) {
        data_[r_offset + c] = static_cast<cost>(src(r, c) * scale_factor);
      }
      
      padRowAfterCol(r, n_col, max_score);
    }
    
    padAfterRow(n_row, max_score);
  }
  
  // With this optimization, the cost of copying the data is almost offset
  // by the savings in reading column-wise.
  void transpose() noexcept {
    const cost* __restrict data_ptr = data_.data();
    cost* __restrict t_data_ptr = t_data_.data();
    
#if defined(__GNUC__) || defined(__clang__)
    data_ptr = static_cast<const cost*>(__builtin_assume_aligned(data_ptr, 64));
    t_data_ptr = static_cast<cost*>(__builtin_assume_aligned(t_data_ptr, 64));
#endif
    for (size_t i = 0; i < dim_; i += BLOCK_SIZE) {
      for (size_t j = 0; j < dim_; j += BLOCK_SIZE) {
        for (size_t r = i; r < std::min(i + BLOCK_SIZE, dim_); ++r) {
          for (size_t c = j; c < std::min(j + BLOCK_SIZE, dim_); ++c) {
            t_data_ptr[c * dim8_ + r] = data_ptr[r * dim8_ + c];
          }
        }
      }
    }
  }
  
  // Access operator for read/write
  cost& operator()(lap_row i, lap_col j) {
    return data_[static_cast<size_t>(i) * dim8_ + j];
  }
  
  // Const version for read-only access
  [[nodiscard]] const cost& operator()(lap_row i, lap_col j) const {
    return data_[static_cast<size_t>(i) * dim8_ + j];
  }
  
  [[nodiscard]] const cost& row0(lap_col j) const {
    return data_[j];
  }
  
  [[nodiscard]] const cost& entry0(lap_row i) const {
    return data_[static_cast<size_t>(i) * dim8_];
  }
  
  [[nodiscard]] cost* row(lap_row i) {
    return &data_[static_cast<size_t>(i) * dim8_];
  }
  
  [[nodiscard]] const cost* row(lap_row i) const {
    return &data_[static_cast<size_t>(i) * dim8_];
  }
  
  [[nodiscard]] cost* col(lap_col i) {
    return &t_data_[static_cast<size_t>(i) * dim8_];
  }
  
  [[nodiscard]] const cost* col(lap_col i) const {
    return &t_data_[static_cast<size_t>(i) * dim8_];
  }
  
  void padAfterRow(lap_row start_row, cost value) {
    size_t start_index = static_cast<size_t>(start_row) * dim8_;
    std::fill(data_.begin() + start_index, data_.end(), value);
  }
  
  void padRowAfterCol(const lap_row r, const lap_col start_col,
                      const cost value) {
    size_t r_offset = r * dim8_;
    size_t actual_start_col = static_cast<size_t>(start_col);
    size_t start_index = r_offset + actual_start_col;
    size_t end_index = start_index + dim_ - actual_start_col;
    std::fill(data_.begin() + start_index, data_.begin() + end_index, value);
  }
  
  std::pair<cost, lap_dim> findColMin(const lap_col j) {
    cost min = firstInCol(j);
    size_t imin = 0;
    const size_t i_limit = (dim_ < 4 ? 0 : dim_ - 3);
    
    for (size_t i = 1; i < i_limit; i += 4) {
      // Very modest performance gain from unrolling this loop.
      const cost c0 = nextInCol();
      const cost c1 = nextInCol();
      const cost c2 = nextInCol();
      const cost c3 = nextInCol();
      
      if (c0 < min) { min = c0; imin = i; }
      if (c1 < min) { min = c1; imin = i + 1; }
      if (c2 < min) { min = c2; imin = i + 2; }
      if (c3 < min) { min = c3; imin = i + 3; }
    }
    
    for (size_t i = i_limit; i < dim_; ++i) {
      const cost current_cost = nextInCol();
      if (current_cost < min) {
        min = current_cost;
        imin = i;
      }
    }
    return {min, static_cast<lap_row>(imin)};
  }
  
};

using cost_matrix = CostMatrix;

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