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
  bool transposed_;
  
  const size_t block_containing(const size_t x) {
    return ((x + BLOCK_SIZE - 1) / BLOCK_SIZE) * BLOCK_SIZE;
  }
  
public:
  CostMatrix(size_t dim)
    : dim_(dim),
      dim8_(block_containing(dim_)),
      data_(std::vector<cost>(dim8_ * dim8_)),
      t_data_(std::vector<cost>(dim8_ * dim8_)),
      transposed_(false) {}
      
  CostMatrix(const Rcpp::NumericMatrix& src, const double x_max)
    : dim_((std::max(src.nrow(), src.ncol()))),  // or pad here as needed
      dim8_(block_containing(dim_)),
      data_(std::vector<cost>(dim8_ * dim_)),
      t_data_(std::vector<cost>(dim8_ * dim_)),
      transposed_(false)
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
  
  void makeTranspose() noexcept {
    if (transposed_) return;
    
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
    
    transposed_ = true;
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
  
  std::pair<cost, lap_row> findColMin(lap_col j) {
    makeTranspose();
    const cost* col_data = col(j);
    const auto min_ptr = std::min_element(col_data, col_data + dim_);
    return {*min_ptr,
            static_cast<lap_row>(std::distance(col_data, min_ptr))};
  }
  
  // test2000 = 207 ms
  std::tuple<cost, cost, lap_col, lap_col> findRowSubmin(
      const lap_row* i, const std::vector<cost>& v
  ) const {
    const cost* __restrict v_ptr = v.data();
    const cost* __restrict row_i = row(*i);
    const lap_col dim = static_cast<lap_col>(dim_);
    
    //     Find minimum and second minimum reduced cost over columns.
    cost min_val = row_i[0] - v_ptr[0];
    lap_col min_idx = 0;
    lap_col submin_idx = 0;
    cost submin_val = BIG;
    
    lap_col j = 1;
    const lap_col j_limit = (dim < 4 ? 0 : static_cast<lap_col>(dim - 3));
    
    for (; j < j_limit; j += 4) {
      assert(BLOCK_SIZE >= 4);  // Unrolling loop x4 gives ~20% speedup
      const cost h0 = row_i[j] - v_ptr[j];
      const cost h1 = row_i[j + 1] - v_ptr[j + 1];
      const cost h2 = row_i[j + 2] - v_ptr[j + 2];
      const cost h3 = row_i[j + 3] - v_ptr[j + 3];
      if (h0 < submin_val) {
        if (h0 < min_val) {
          submin_val = min_val;
          min_val = h0;
          submin_idx = min_idx;
          min_idx = j;
        } else {
          submin_val = h0;
          submin_idx = j;
        }
      }
      if (h1 < submin_val) {
        if (h1 < min_val) {
          submin_val = min_val;
          min_val = h1;
          submin_idx = min_idx;
          min_idx = j + 1;
        } else {
          submin_val = h1;
          submin_idx = j + 1;
        }
      }
      if (h2 < submin_val) {
        if (h2 < min_val) {
          submin_val = min_val;
          min_val = h2;
          submin_idx = min_idx;
          min_idx = j + 2;
        } else {
          submin_val = h2;
          submin_idx = j + 2;
        }
      }
      if (h3 < submin_val) {
        if (h3 < min_val) {
          submin_val = min_val;
          min_val = h3;
          submin_idx = min_idx;
          min_idx = j + 3;
        } else {
          submin_val = h3;
          submin_idx = j + 3;
        }
      }
    }
    for (; j < dim; ++j) {
      const cost h = row_i[j] - v_ptr[j];
      if (h < submin_val) {
        if (h < min_val) {
          submin_val = min_val;
          min_val = h;
          submin_idx = min_idx;
          min_idx = j;
        } else {
          submin_val = h;
          submin_idx = j;
        }
      }
    }
    return {min_val, submin_val, min_idx, submin_idx};
  }
  
  // test2000 = 260 ms
  std::tuple<cost, cost, lap_col, lap_col> findRowSubminNaive(
      const lap_row* i,
      const std::vector<cost>& v_ptr
  ) {
    const cost* __restrict row_i = row(*i);
    const lap_col dim = static_cast<lap_col>(dim_);
    
    if (dim < 2) {
      return {row_i[0] - v_ptr[0], std::numeric_limits<cost>::max(), 0, 0};
    }
    
    // Initialize with first two elements
    cost h0 = row_i[0] - v_ptr[0];
    cost h1 = row_i[1] - v_ptr[1];
    
    cost min_val, submin_val;
    lap_col min_idx, submin_idx;
    
    if (h0 <= h1) {
      min_val = h0; submin_val = h1;
      min_idx = 0; submin_idx = 1;
    } else {
      min_val = h1; submin_val = h0;
      min_idx = 1; submin_idx = 0;
    }
    
    // Process remaining elements
    for (lap_col j = 2; j < dim; ++j) {
      const cost h = row_i[j] - v_ptr[j];
      
      if (h < min_val) {
        submin_val = min_val;
        submin_idx = min_idx;
        min_val = h;
        min_idx = j;
      } else if (h < submin_val) {
        submin_val = h;
        submin_idx = j;
      }
    }
    
    return {min_val, submin_val, min_idx, submin_idx};
  }
  
  // test2000 = 370 ms (!)
  std::tuple<cost, cost, lap_col, lap_col> findRowSubminTwoPassNaive(
      const lap_row* i,
      const std::vector<cost>& v_ptr
  ) {
    const cost* __restrict row_i = row(*i);
    cost min_val = std::numeric_limits<cost>::max();
    lap_col min_idx = 0;
    const lap_col dim = static_cast<lap_col>(dim_);
    
    for (lap_col j = 0; j < dim; ++j) {
      const cost h = row_i[j] - v_ptr[j];
      if (h < min_val) {
        min_val = h;
        min_idx = j;
      }
    }
    
    // Second pass: find subminimum
    cost submin_val = std::numeric_limits<cost>::max();
    lap_col submin_idx = (min_idx == 0) ? 1 : 0;
    
    for (lap_col j = 0; j < dim; ++j) {
      if (j != min_idx) {
        const cost h = row_i[j] - v_ptr[j];
        if (h < submin_val) {
          submin_val = h;
          submin_idx = j;
        }
      }
    }
    
    return {min_val, submin_val, min_idx, submin_idx};
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