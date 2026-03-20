#ifndef TREEDIST_LAP_IMPL_H_
#define TREEDIST_LAP_IMPL_H_

// LAP (Linear Assignment Problem) — Jonker–Volgenant implementation.
//
// Guard this with #define TREEDIST_LAP_IMPLEMENTATION before including.
// Include in exactly one translation unit per package.
//
// Interrupt handling:  define TREEDIST_CHECK_INTERRUPT() before including
// to enable user-interrupt checks (e.g., Rcpp::checkUserInterrupt()).
// If not defined, interrupt checks are silently skipped.
//
// Original algorithm:
//   R. Jonker and A. Volgenant, "A Shortest Augmenting Path Algorithm
//   for Dense and Sparse Linear Assignment Problems," Computing 38,
//   325-340, 1987.

#ifdef TREEDIST_LAP_IMPLEMENTATION

#include "lap.h"

#ifndef TREEDIST_CHECK_INTERRUPT
#define TREEDIST_CHECK_INTERRUPT() ((void)0)
#define TREEDIST_CHECK_INTERRUPT_DEFINED_HERE_
#endif

namespace TreeDist {

namespace detail {
  inline bool nontrivially_less_than(cost a, cost b) noexcept {
    return a + ((a > ROUND_PRECISION) ? 8 : 0) < b;
  }
} // namespace detail

cost lap(const lap_row dim,
         CostMatrix& input_cost,
         std::vector<lap_col>& rowsol,
         std::vector<lap_row>& colsol,
         const bool allow_interrupt,
         LapScratch& scratch)
{
  lap_row num_free = 0;
  scratch.ensure(dim);
  auto& v       = scratch.v;
  auto& matches = scratch.matches;
  std::fill(matches.begin(), matches.begin() + dim, 0);
  const cost* __restrict__ v_ptr = v.data();

  // COLUMN REDUCTION
  for (lap_col j = dim; j--; ) {
    const auto [min, imin] = input_cost.findColMin(j);
    v[j] = min;
    ++matches[imin];

    if (matches[imin] == 1) {
      rowsol[imin] = j;
      colsol[j] = imin;
    } else if (v_ptr[j] < v_ptr[rowsol[imin]]) {
      const lap_col j1 = rowsol[imin];
      rowsol[imin] = j;
      colsol[j] = imin;
      colsol[j1] = -1;
    } else {
      colsol[j] = -1;
    }
  }

  // REDUCTION TRANSFER
  auto& freeunassigned = scratch.freeunassigned;

  for (lap_row i = 0; i < dim; ++i) {
    if (matches[i] == 0) {
      freeunassigned[num_free++] = i;
    } else if (matches[i] == 1) {
      const lap_col j1 = rowsol[i];
      const cost* row_i = input_cost.row(i);
      cost min_cost;
      if (j1 == 0) {
        min_cost = row_i[1] - v_ptr[1];
        for (lap_col j = 2; j < dim; ++j) {
          const cost reduced_cost = row_i[j] - v_ptr[j];
          if (reduced_cost < min_cost) {
            min_cost = reduced_cost;
          }
        }
      } else {
        min_cost = row_i[0] - v_ptr[0];
        for (lap_col j = 1; j < dim; ++j) {
          if (j == j1) continue;
          const cost reduced_cost = row_i[j] - v_ptr[j];
          if (reduced_cost < min_cost) {
            min_cost = reduced_cost;
          }
        }
      }
      v[j1] -= min_cost;
    }
  }

  // AUGMENTING ROW REDUCTION
  auto& col_list = scratch.col_list;
  int loopcnt = 0;

  do {
    ++loopcnt;

    lap_row previous_num_free = num_free;
    num_free = 0;
    lap_row k = 0;
    while (k < previous_num_free) {
      const lap_row i = freeunassigned[k++];
      const auto [umin, usubmin, min_idx, j2] =
        input_cost.findRowSubmin(&i, v);
      lap_col j1 = min_idx;

      lap_row i0 = colsol[j1];
      const bool strictly_less =
        detail::nontrivially_less_than(umin, usubmin);
      if (strictly_less) {
        v[j1] -= (usubmin - umin);
      } else if (i0 > -1) {
        j1 = j2;
        i0 = colsol[j2];
      }

      rowsol[i] = j1;
      colsol[j1] = i;

      if (i0 > -1) {
        if (strictly_less) {
          freeunassigned[--k] = i0;
          if (allow_interrupt) TREEDIST_CHECK_INTERRUPT();
        } else {
          freeunassigned[num_free++] = i0;
        }
      }
    }
  } while (loopcnt < 2);

  // AUGMENT SOLUTION for each free row.
  auto& d           = scratch.d;
  auto& predecessor = scratch.predecessor;

  for (lap_row f = 0; f < num_free; ++f) {
    bool unassignedfound = false;
    lap_row free_row = freeunassigned[f];
    const cost* free_row_cost = input_cost.row(free_row);
    lap_col endofpath = 0;
    lap_col last = 0;
    lap_row i;
    lap_col j1;

    for (lap_col j = 0; j < dim; ++j) {
      d[j] = free_row_cost[j] - v_ptr[j];
      predecessor[j] = free_row;
      col_list[j] = j;
    }

    cost min = 0;
    lap_col low = 0;
    lap_col up = 0;

    do {
      if (up == low) {
        last = low - 1;
        min = d[col_list[up++]];

        for (lap_dim k = up; k < dim; ++k) {
          const lap_col j = col_list[k];
          const cost h = d[j];
          if (h <= min) {
            if (h < min) {
              up = low;
              min = h;
            }
            col_list[k] = col_list[up];
            col_list[up++] = j;
          }
        }
        for (lap_dim k = low; k < up; ++k) {
          if (colsol[col_list[k]] < 0) {
            endofpath = col_list[k];
            unassignedfound = true;
            break;
          }
        }
      }

      if (!unassignedfound) {
        j1 = col_list[low++];
        i = colsol[j1];
        const cost* row_i = input_cost.row(i);
        const cost h = row_i[j1] - v_ptr[j1] - min;

        for (lap_dim k = up; k < dim; ++k) {
          const lap_col j = col_list[k];
          cost v2 = row_i[j] - v_ptr[j] - h;
          if (v2 < d[j]) {
            predecessor[j] = i;
            if (v2 == min) {
              if (colsol[j] < 0) {
                endofpath = j;
                unassignedfound = true;
                break;
              } else {
                col_list[k] = col_list[up];
                col_list[up++] = j;
              }
            }
            d[j] = v2;
          }
        }
      }
    } while (!unassignedfound);

    for (lap_dim k = 0; k <= last; ++k) {
      j1 = col_list[k];
      v[j1] += d[j1] - min;
    }

    do {
      i = predecessor[endofpath];
      colsol[endofpath] = i;
      j1 = endofpath;
      endofpath = rowsol[i];
      rowsol[i] = j1;
    } while (i != free_row);
  }

  // Calculate optimal cost.
  cost lapcost = 0;
  for (lap_dim i = 0; i < dim; ++i) {
    lapcost += input_cost(i, rowsol[i]);
  }

  return lapcost;
}

} // namespace TreeDist

#ifdef TREEDIST_CHECK_INTERRUPT_DEFINED_HERE_
#undef TREEDIST_CHECK_INTERRUPT
#undef TREEDIST_CHECK_INTERRUPT_DEFINED_HERE_
#endif

#endif // TREEDIST_LAP_IMPLEMENTATION
#endif // TREEDIST_LAP_IMPL_H_
