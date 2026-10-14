// Copy of TreeDist::lap() (inst/include/TreeDist/lap_impl.h) with
//  - a hook called after each phase (with the current column duals v), and
//  - an optional warm start from given column duals.
// Phases: 0 = column reduction, 1 = reduction transfer, 2 = ARR,
//         3 = after each augmentation (arg = augmentations done, of total).
#pragma once
#include "lap.h"

namespace ana {
using namespace TreeDist;

inline bool nontrivially_less_than(cost a, cost b) noexcept {
  return a + ((a > ROUND_PRECISION) ? 8 : 0) < b;
}

struct LapWork {
  long long n_free_after_init = 0; // free rows entering ARR
  long long n_free_after_arr = 0;  // free rows entering augmentation
  long long arr_steps = 0;         // findRowSubmin calls
  long long aug_inner = 0;         // Dijkstra inner-loop column visits
  long long aug_scan = 0;          // Dijkstra min-scan column visits
  unsigned long long t_init = 0, t_arr = 0, t_aug = 0;
};

template <typename Hook>
cost lap_hooked(const lap_row dim, CostMatrix &input_cost,
                std::vector<lap_col> &rowsol, std::vector<lap_row> &colsol,
                LapScratch &scratch, const cost *v_init, Hook &&hook,
                LapWork &work, const int arr_loops = 2) {
  unsigned long long tick = __rdtsc();
  lap_row num_free = 0;
  scratch.ensure(dim);
  auto &v = scratch.v;
  auto &matches = scratch.matches;
  std::fill(matches.begin(), matches.begin() + dim, 0);
  const cost *__restrict__ v_ptr = v.data();
  auto &freeunassigned = scratch.freeunassigned;

  if (!v_init) {
    // COLUMN REDUCTION
    for (lap_col j = dim; j--;) {
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
    hook(0, v, 0, 0);

    // REDUCTION TRANSFER
    for (lap_row i = 0; i < dim; ++i) {
      if (matches[i] == 0) {
        freeunassigned[num_free++] = i;
      } else if (matches[i] == 1) {
        const lap_col j1 = rowsol[i];
        const cost *row_i = input_cost.row(i);
        cost min_cost;
        if (j1 == 0) {
          min_cost = row_i[1] - v_ptr[1];
          for (lap_col j = 2; j < dim; ++j) {
            const cost reduced_cost = row_i[j] - v_ptr[j];
            if (reduced_cost < min_cost) min_cost = reduced_cost;
          }
        } else {
          min_cost = row_i[0] - v_ptr[0];
          for (lap_col j = 1; j < dim; ++j) {
            if (j == j1) continue;
            const cost reduced_cost = row_i[j] - v_ptr[j];
            if (reduced_cost < min_cost) min_cost = reduced_cost;
          }
        }
        v[j1] -= min_cost;
      }
    }
    hook(1, v, 0, 0);
  } else {
    // WARM START: v given; assign each row to its cheapest reduced column
    // if that column is still free, so assigned edges are tight.
    for (lap_col j = 0; j < dim; ++j) {
      v[j] = v_init[j];
      colsol[j] = -1;
    }
    for (lap_row i = 0; i < dim; ++i) {
      const cost *row_i = input_cost.row(i);
      cost best = row_i[0] - v_ptr[0];
      lap_col best_j = 0;
      for (lap_col j = 1; j < dim; ++j) {
        const cost h = row_i[j] - v_ptr[j];
        if (h < best) { best = h; best_j = j; }
      }
      if (colsol[best_j] < 0) {
        colsol[best_j] = i;
        rowsol[i] = best_j;
      } else {
        freeunassigned[num_free++] = i;
      }
    }
    hook(1, v, 0, 0);
  }
  work.n_free_after_init += num_free;
  { unsigned long long t = __rdtsc(); work.t_init += t - tick; tick = t; }

  // AUGMENTING ROW REDUCTION
  int loopcnt = 0;
  while (loopcnt < arr_loops) {
    ++loopcnt;
    lap_row previous_num_free = num_free;
    num_free = 0;
    lap_row k = 0;
    while (k < previous_num_free) {
      const lap_row i = freeunassigned[k++];
      ++work.arr_steps;
      const auto [umin, usubmin, min_idx, j2] = input_cost.findRowSubmin(&i, v);
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
  }
  hook(2, v, 0, num_free);
  work.n_free_after_arr += num_free;
  { unsigned long long t = __rdtsc(); work.t_arr += t - tick; tick = t; }

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
          ++work.aug_scan;
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
          ++work.aug_inner;
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
    hook(3, v, f + 1, num_free);
  }

  { unsigned long long t = __rdtsc(); work.t_aug += t - tick; tick = t; }
  cost lapcost = 0;
  for (lap_dim i = 0; i < dim; ++i) {
    lapcost += input_cost(i, rowsol[i]);
  }
  return lapcost;
}
} // namespace ana
