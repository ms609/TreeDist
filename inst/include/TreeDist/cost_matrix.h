#ifndef TREEDIST_COST_MATRIX_H_
#define TREEDIST_COST_MATRIX_H_

// TreeDist CostMatrix — cache-aligned, block-padded cost matrix for the LAP.
//
// Rcpp-free public API.  The Rcpp-dependent constructor
// (CostMatrix(Rcpp::NumericMatrix)) lives in TreeDist's own src/lap.h.

#include "types.h"

#include <algorithm>
#include <cassert>
#include <tuple>
#include <vector>

namespace TreeDist {

class CostMatrix {
private:
  std::size_t dim_;
  std::size_t dim8_;
  alignas(64) std::vector<cost> data_;
  alignas(64) std::vector<cost> t_data_;
  bool transposed_;

  static std::size_t block_containing(const std::size_t x) {
    return ((x + BLOCK_SIZE - 1) / BLOCK_SIZE) * BLOCK_SIZE;
  }

public:
  // Default constructor for pooled instances (zero-size, no allocation).
  CostMatrix() : dim_(0), dim8_(0), transposed_(false) {}

  explicit CostMatrix(std::size_t dim)
    : dim_(dim),
      dim8_(block_containing(dim_)),
      data_(std::vector<cost>(dim8_ * dim8_)),
      t_data_(std::vector<cost>(dim8_ * dim8_)),
      transposed_(false) {}

  // Resize for reuse.  Only reallocates when the new dimension exceeds the
  // current buffer capacity; otherwise just updates dim_/dim8_ and marks
  // the transpose as stale.
  void resize(std::size_t new_dim) {
    dim_ = new_dim;
    dim8_ = block_containing(new_dim);
    const std::size_t needed = dim8_ * dim8_;
    if (data_.size() < needed) {
      data_.resize(needed);
      t_data_.resize(needed);
    }
    transposed_ = false;
  }

  void reset() noexcept { transposed_ = false; }

  [[nodiscard]] std::size_t dim() const noexcept { return dim_; }

  // ---- Element access ----

  cost& operator()(lap_row i, lap_col j) {
    return data_[static_cast<std::size_t>(i) * dim8_ + j];
  }

  [[nodiscard]] const cost& operator()(lap_row i, lap_col j) const {
    return data_[static_cast<std::size_t>(i) * dim8_ + j];
  }

  [[nodiscard]] const cost& row0(lap_col j) const {
    return data_[j];
  }

  [[nodiscard]] const cost& entry0(lap_row i) const {
    return data_[static_cast<std::size_t>(i) * dim8_];
  }

  [[nodiscard]] cost* row(lap_row i) {
    return &data_[static_cast<std::size_t>(i) * dim8_];
  }

  [[nodiscard]] const cost* row(lap_row i) const {
    return &data_[static_cast<std::size_t>(i) * dim8_];
  }

  [[nodiscard]] cost* col(lap_col i) {
    return &t_data_[static_cast<std::size_t>(i) * dim8_];
  }

  [[nodiscard]] const cost* col(lap_col i) const {
    return &t_data_[static_cast<std::size_t>(i) * dim8_];
  }

  // Write a value to both the data and transposed-data arrays.
  // After filling the entire matrix this way, call markTransposed()
  // to skip the lazy transpose in makeTranspose().
  void setWithTranspose(lap_row i, lap_col j, cost value) {
    data_[static_cast<std::size_t>(i) * dim8_ + j] = value;
    t_data_[static_cast<std::size_t>(j) * dim8_ + i] = value;
  }

  void markTransposed() noexcept { transposed_ = true; }

  // ---- Transpose ----

  void makeTranspose() noexcept {
    if (transposed_) return;

    const cost* __restrict__ data_ptr = data_.data();
    cost* __restrict__ t_data_ptr = t_data_.data();

#if defined(__GNUC__) || defined(__clang__)
    data_ptr =
      static_cast<const cost*>(__builtin_assume_aligned(data_ptr, 64));
    t_data_ptr =
      static_cast<cost*>(__builtin_assume_aligned(t_data_ptr, 64));
#endif
    for (std::size_t i = 0; i < dim_; i += BLOCK_SIZE) {
      for (std::size_t j = 0; j < dim_; j += BLOCK_SIZE) {
        for (std::size_t r = i;
             r < std::min(i + BLOCK_SIZE, dim_); ++r) {
          for (std::size_t c = j;
               c < std::min(j + BLOCK_SIZE, dim_); ++c) {
            t_data_ptr[c * dim8_ + r] = data_ptr[r * dim8_ + c];
          }
        }
      }
    }
    transposed_ = true;
  }

  void makeUntranspose() noexcept {
    const cost* __restrict__ data_ptr = t_data_.data();
    cost* __restrict__ t_data_ptr = data_.data();

#if defined(__GNUC__) || defined(__clang__)
    data_ptr =
      static_cast<const cost*>(__builtin_assume_aligned(data_ptr, 64));
    t_data_ptr =
      static_cast<cost*>(__builtin_assume_aligned(t_data_ptr, 64));
#endif
    for (std::size_t i = 0; i < dim_; i += BLOCK_SIZE) {
      for (std::size_t j = 0; j < dim_; j += BLOCK_SIZE) {
        for (std::size_t r = i;
             r < std::min(i + BLOCK_SIZE, dim_); ++r) {
          for (std::size_t c = j;
               c < std::min(j + BLOCK_SIZE, dim_); ++c) {
            t_data_ptr[c * dim8_ + r] = data_ptr[r * dim8_ + c];
          }
        }
      }
    }
  }

  // ---- Padding ----

  void padTrAfterCol(lap_row start_row, cost value) {
    std::size_t start_index = static_cast<std::size_t>(start_row) * dim8_;
    std::size_t end_index   = dim_ * dim8_;
    std::fill(t_data_.begin() + start_index,
              t_data_.begin() + end_index, value);
  }

  void padAfterRow(lap_row start_row, cost value) {
    std::size_t start_index = static_cast<std::size_t>(start_row) * dim8_;
    std::size_t end_index   = dim_ * dim8_;
    std::fill(data_.begin() + start_index,
              data_.begin() + end_index, value);
  }

  void padTrColAfterRow(const lap_row r, const lap_col start_col,
                        const cost value) {
    std::size_t r_offset = r * dim8_;
    std::size_t actual_start_col = static_cast<std::size_t>(start_col);
    std::size_t start_index = r_offset + actual_start_col;
    std::size_t end_index = start_index + dim_ - actual_start_col;
    std::fill(t_data_.begin() + start_index,
              t_data_.begin() + end_index, value);
  }

  void padRowAfterCol(const lap_row r, const lap_col start_col,
                      const cost value) {
    std::size_t r_offset = r * dim8_;
    std::size_t actual_start_col = static_cast<std::size_t>(start_col);
    std::size_t start_index = r_offset + actual_start_col;
    std::size_t end_index = start_index + dim_ - actual_start_col;
    std::fill(data_.begin() + start_index,
              data_.begin() + end_index, value);
  }

  // ---- Search ----

  std::pair<cost, lap_row> findColMin(lap_col j,
                                      lap_row search_dim = -1) {
    makeTranspose();
    const cost* col_data = col(j);
    const lap_row n =
      (search_dim < 0) ? static_cast<lap_row>(dim_) : search_dim;
    const auto min_ptr = std::min_element(col_data, col_data + n);
    return {*min_ptr,
            static_cast<lap_row>(std::distance(col_data, min_ptr))};
  }

  std::tuple<cost, cost, lap_col, lap_col> findRowSubmin(
      const lap_row* i, const std::vector<cost>& v
  ) const {
    assert(dim_ > 1);

    const cost* __restrict__ row_i = row(*i);
    const lap_col dim = static_cast<lap_col>(dim_);
    const cost* __restrict__ v_ptr = v.data();

    const cost h0 = row_i[0] - v_ptr[0];
    const cost h1 = row_i[1] - v_ptr[1];

    cost min_val, submin_val;
    lap_col min_idx, submin_idx;

    if (h0 > h1) {
      min_val = h1; submin_val = h0;
      min_idx = 1;  submin_idx = 0;
    } else {
      min_val = h0; submin_val = h1;
      min_idx = 0;  submin_idx = 1;
    }

    const lap_col j_limit =
      (dim < 4 ? 0 : static_cast<lap_col>(dim - 3));

    for (lap_col j = 2; j < j_limit; j += 4) {
      assert(BLOCK_SIZE >= 4);
      const cost h0 = row_i[j]     - v_ptr[j];
      const cost h1 = row_i[j + 1] - v_ptr[j + 1];
      const cost h2 = row_i[j + 2] - v_ptr[j + 2];
      const cost h3 = row_i[j + 3] - v_ptr[j + 3];
      if (h0 < submin_val) {
        if (h0 < min_val) {
          submin_val = min_val; min_val = h0;
          submin_idx = min_idx; min_idx = j;
        } else { submin_val = h0; submin_idx = j; }
      }
      if (h1 < submin_val) {
        if (h1 < min_val) {
          submin_val = min_val; min_val = h1;
          submin_idx = min_idx; min_idx = j + 1;
        } else { submin_val = h1; submin_idx = j + 1; }
      }
      if (h2 < submin_val) {
        if (h2 < min_val) {
          submin_val = min_val; min_val = h2;
          submin_idx = min_idx; min_idx = j + 2;
        } else { submin_val = h2; submin_idx = j + 2; }
      }
      if (h3 < submin_val) {
        if (h3 < min_val) {
          submin_val = min_val; min_val = h3;
          submin_idx = min_idx; min_idx = j + 3;
        } else { submin_val = h3; submin_idx = j + 3; }
      }
    }
    for (lap_col j = j_limit; j < dim; ++j) {
      const cost h = row_i[j] - v_ptr[j];
      if (h < submin_val) {
        if (h < min_val) {
          submin_val = min_val; min_val = h;
          submin_idx = min_idx; min_idx = j;
        } else { submin_val = h; submin_idx = j; }
      }
    }
    return {min_val, submin_val, min_idx, submin_idx};
  }
};

} // namespace TreeDist

#endif // TREEDIST_COST_MATRIX_H_
