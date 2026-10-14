// Jonker-Volgenant LAP warm-started from given column duals.
//
// As TreeDist::lap() (inst/include/TreeDist/lap_impl.h), except that column
// reduction and reduction transfer are replaced by a start from the column
// duals `v_init`: each row takes its cheapest reduced-cost column if that
// column is still free, so every assigned edge is tight, as the augmenting
// row reduction and augmentation phases require.  Any `v_init` gives the
// same optimum; duals close to the optimum leave few rows to augment.
//
// Kept in its own translation unit so that lap()'s code generation, which is
// sensitive to layout, is unaffected.  Final column duals are left in
// scratch.v.

#include "lap.h"
#include <algorithm>

namespace {
inline bool nontrivially_less_than(cost a, cost b) noexcept {
  return a + ((a > ROUND_PRECISION) ? 8 : 0) < b;
}
}

cost lap_warm(const lap_row dim, CostMatrix &input_cost,
              std::vector<lap_col> &rowsol, std::vector<lap_row> &colsol,
              LapScratch &scratch, const cost *v_init) {
  lap_row num_free = 0;
  scratch.ensure(dim);
  auto &v = scratch.v;
  auto &freeunassigned = scratch.freeunassigned;
  const cost *__restrict__ v_ptr = v.data();

  // WARM START
  for (lap_col j = 0; j < dim; ++j) {
    v[j] = v_init[j];
    colsol[j] = -1;
  }
  for (lap_row i = 0; i < dim; ++i) {
    const cost *__restrict__ row_i = input_cost.row(i);
    cost best = row_i[0] - v_ptr[0];
    lap_col best_j = 0;
    for (lap_col j = 1; j < dim; ++j) {
      const cost h = row_i[j] - v_ptr[j];
      if (h < best) {
        best = h;
        best_j = j;
      }
    }
    if (colsol[best_j] < 0) {
      colsol[best_j] = i;
      rowsol[i] = best_j;
    } else {
      freeunassigned[num_free++] = i;
    }
  }

  // AUGMENTING ROW REDUCTION
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
      const bool strictly_less = nontrivially_less_than(umin, usubmin);
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
        } else {
          freeunassigned[num_free++] = i0;
        }
      }
    }
  } while (loopcnt < 2);

  // AUGMENT SOLUTION for each free row.
  cost *__restrict__ d_ptr = scratch.d.data();
  lap_row *__restrict__ pred_ptr = scratch.predecessor.data();
  lap_col *__restrict__ cl_ptr = scratch.col_list.data();

  for (lap_row f = 0; f < num_free; ++f) {
    bool unassignedfound = false;
    lap_row free_row = freeunassigned[f];
    const cost *__restrict__ free_row_cost = input_cost.row(free_row);
    lap_col endofpath = 0;
    lap_col last = 0;
    lap_row i;
    lap_col j1;

    for (lap_col j = 0; j < dim; ++j) {
      d_ptr[j] = free_row_cost[j] - v_ptr[j];
      pred_ptr[j] = free_row;
      cl_ptr[j] = j;
    }

    cost min = 0;
    lap_col low = 0;
    lap_col up = 0;

    do {
      if (up == low) {
        last = low - 1;
        min = d_ptr[cl_ptr[up++]];
        for (lap_dim k = up; k < dim; ++k) {
          const lap_col j = cl_ptr[k];
          const cost h = d_ptr[j];
          if (h <= min) {
            if (h < min) {
              up = low;
              min = h;
            }
            cl_ptr[k] = cl_ptr[up];
            cl_ptr[up++] = j;
          }
        }
        for (lap_dim k = low; k < up; ++k) {
          if (colsol[cl_ptr[k]] < 0) {
            endofpath = cl_ptr[k];
            unassignedfound = true;
            break;
          }
        }
      }

      if (!unassignedfound) {
        j1 = cl_ptr[low++];
        i = colsol[j1];
        const cost *__restrict__ row_i = input_cost.row(i);
        const cost h = row_i[j1] - v_ptr[j1] - min;
        for (lap_dim k = up; k < dim; ++k) {
          const lap_col j = cl_ptr[k];
          cost v2 = row_i[j] - v_ptr[j] - h;
          if (v2 < d_ptr[j]) {
            pred_ptr[j] = i;
            if (v2 == min) {
              if (colsol[j] < 0) {
                endofpath = j;
                unassignedfound = true;
                break;
              } else {
                cl_ptr[k] = cl_ptr[up];
                cl_ptr[up++] = j;
              }
            }
            d_ptr[j] = v2;
          }
        }
      }
    } while (!unassignedfound);

    for (lap_dim k = 0; k <= last; ++k) {
      j1 = cl_ptr[k];
      v[j1] += d_ptr[j1] - min;
    }

    do {
      i = pred_ptr[endofpath];
      colsol[endofpath] = i;
      j1 = endofpath;
      endofpath = rowsol[i];
      rowsol[i] = j1;
    } while (i != free_row);
  }

  cost lapcost = 0;
  for (lap_dim i = 0; i < dim; ++i) {
    lapcost += input_cost(i, rowsol[i]);
  }
  return lapcost;
}
