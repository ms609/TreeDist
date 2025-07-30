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
 CHANGED 2020-01-01 by Martin Smith <martin.smith@durham.ac.uk> for 
 integration with R.
 CHANGED 2016-05-13 by Yong Yang(yongyanglink@gmail.com) in column reduction 
 part according to matlab version of LAPJV algorithm 
 https://github.com/yongyanghz/LAPJV-algorithm-c/blob/master/LAPJV/lap.cpp
 (Copyright (c) 2010, Yi Cao All rights reserved)--
 https://www.mathworks.com/matlabcentral/fileexchange/26836-lapjv-jonker-volgenant-algorithm-for-linear-assignment-problem-v3-0:
 *
 *************************************************************************/
#include "tree_distances.h"
using namespace Rcpp;


// [[Rcpp::export]]
List lapjv(NumericMatrix x, NumericVector maxX) {
  const int16 n_row = x.nrow();
  const int16 n_col = x.ncol();
  const int16 max_dim = (n_row > n_col) ? n_row : n_col;
  const int16 spare_rows = n_row - n_col;
  const cost max_score = cost(BIG / max_dim);
  
  std::vector<lap_col> rowsol(max_dim);
  std::vector<lap_row> colsol(max_dim);
  std::vector<cost> u(max_dim);
  std::vector<cost> v(max_dim);
  
  cost_matrix input(max_dim);
  
  const double x_max = maxX[0];
  const double scale_factor = max_score / x_max;
  for (int16 r = 0; r < n_row; ++r) {
    const row_offset r_offset = static_cast<size_t>(r) * max_dim;
    for (int16 c = 0; c < n_col; ++c) {
      input.fromRow(r_offset, c) = cost(x(r, c) * scale_factor);
    }
    input.padRowAfterCol(r_offset, n_col, max_score);
    for (int16 c = n_col; c < max_dim; ++c) {
      input.fromRow(r_offset, c) = max_score;
    }
  }
  input.padAfterRow(n_row, max_score);
  
  cost score = lap(max_dim, input, rowsol, colsol, u, v);
  IntegerVector matching(n_row);
  for (int16 i = 0; i < n_row; ++i) {
    matching[i] = rowsol[i] < n_col ? rowsol[i] + 1 : NA_INTEGER;
  }
  
  return List::create(
    Named("score") = (double(score) - (std::abs(spare_rows) * max_score))
    / max_score * x_max,
    _["matching"] = matching
  );
}

inline bool nontrivially_less_than(cost a, cost b) noexcept {
  return a + ((a > ROUND_PRECISION) ? 8 : 0) < b;
}

/* This function is the jv shortest augmenting path algorithm to solve the 
   assignment problem */
cost lap(const lap_row dim,
         cost_matrix &input_cost,
         std::vector<lap_col> &rowsol,
         std::vector<lap_row> &colsol,
         std::vector<cost> &u,
         std::vector<cost> &v)
  
  // input:
  // dim        - problem size
  // input_cost - cost matrix
  
  // output:
  // rowsol     - column assigned to row in solution
  // colsol     - row assigned to column in solution
  // u          - dual variables, row reduction numbers
  // v          - dual variables, column reduction numbers
  
{
  bool unassignedfound;
  lap_dim i;
  lap_row num_free = 0;
  
  lap_col j1;
  
  std::vector<lap_col> matches(dim);     // Counts how many times a row could be assigned.
  
  // Init how many times a row will be assigned in the column reduction.
  for (lap_row i = 0; i != dim; i++) {
    
  }
  
  // COLUMN REDUCTION
  for (lap_col j = dim; j--; ) { // Reverse order gives better results.
    // Find minimum cost over rows.
    cost min = input_cost.row0(j);
    lap_dim imin = 0;
    for (lap_row i = 1; i < dim; ++i) {
      const cost current_cost = input_cost(i, j);
      if (current_cost < min) {
        min = current_cost;
        imin = i;
      }
    }
    
    v[j] = min;
    ++matches[imin];
    
    if (matches[imin] == 1) {
      // Init assignment if minimum row assigned for first time.
      rowsol[imin] = j;
      colsol[j] = imin;
    } else if (v[j] < v[rowsol[imin]]) {
      const lap_col j1 = rowsol[imin];
      rowsol[imin] = j;
      colsol[j] = imin;
      colsol[j1] = -1;
    } else {
      colsol[j] = -1; // Row already assigned, column not assigned.
    }
  }
  
  // REDUCTION TRANSFER
  std::vector<lap_row> freeunassigned(dim);        // List of unassigned rows.
  
  for (lap_row i = 0; i < dim; ++i) {
    cost min = 0;
    if (matches[i] == 0) { // Fill list of unassigned 'free' rows.
      freeunassigned[num_free++] = i;
    } else if (matches[i] == 1) { // Transfer reduction from rows that are assigned once.
      const lap_col j1 = rowsol[i];
      min = BIG;
      const row_offset i_offset = static_cast<size_t>(i) * dim;
      for (lap_col j = 0; j < dim; ++j) {
        if (j != j1) {
          const cost reduced_cost = input_cost.fromRow(i_offset, j) - v[j];
          if (reduced_cost < min) {
            min = reduced_cost;
          }
        }
      }
      v[j1] -= min;
    }
  }
  
  //   AUGMENTING ROW REDUCTION
  std::vector<lap_col> col_list(dim);    // List of columns to be scanned in various ways.
  int16 loopcnt = 0;                     // do-loop to be done twice.
  
  
  do {
    ++loopcnt;
    
    //     Scan all free rows.
    //     In some cases, a free row may be replaced with another one to be 
    //     scanned next.
    cost usubmin;
    cost umin;
    lap_row k = 0;
    lap_row previous_num_free = num_free;
    num_free = 0;             // Start list of rows still free after augmenting
                              // row reduction.
    while (k < previous_num_free) {
      const lap_row iii = freeunassigned[k++];
      const row_offset i_offset = static_cast<size_t>(iii) * dim;
      lap_col j2 = 0;
      
      //     Find minimum and second minimum reduced cost over columns.
      umin = input_cost.entry0(static_cast<size_t>(iii)) - v[0];
      j1 = 0;
      usubmin = BIG;
      
      for (lap_col j = 1; j < dim; ++j) {
        const cost h = input_cost(iii, j) - v[j];
        if (h < usubmin) {
          if (h >= umin) {
            usubmin = h;
            j2 = j;
          } else {
            usubmin = umin;
            umin = h;
            j2 = j1;
            j1 = j;
          }
        }
      }
      
      lap_row i0 = colsol[j1];
      if (nontrivially_less_than(umin, usubmin)) {
        //  Change the reduction of the minimum column to increase the minimum
        //  reduced cost in the row to the subminimum.
        v[j1] -= (usubmin - umin);
      } else if (i0 > -1) {
          // Minimum and subminimum equal.
          // Minimum column j1 is assigned.
          // Swap columns j1 and j2, as j2 may be unassigned.
          j1 = j2;
          i0 = colsol[j2];
      }
        
      //    (Re-)assign i to j1, possibly de-assigning an i0.
      rowsol[iii] = j1;
      colsol[j1] = iii;
      
      if (i0 > -1) { // Minimum column j1 assigned earlier.
        if (nontrivially_less_than(umin, usubmin)) {
          // Put in current k, and go back to that k.
          // Continue augmenting path i - j1 with i0.
          freeunassigned[--k] = i0;
          Rcpp::checkUserInterrupt();
        } else {
          // No further augmenting reduction possible.
          // Store i0 in list of free rows for next phase.
          freeunassigned[num_free++] = i0;
        }
      }
    }
  } while (loopcnt < 2); // Repeat once.
  
  // AUGMENT SOLUTION for each free row.
  std::vector<cost> d(dim);              // 'Cost-distance' in augmenting path calculation.
  std::vector<lap_row> predecessor(dim); // Row-predecessor of column in augmenting/alternating path.
  
  for (lap_row f = 0; f < num_free; ++f) {
    lap_row free_row = freeunassigned[f];       // Start row of augmenting path.
    lap_col endofpath = 0;
    lap_col last = 0;
    
    // Dijkstra shortest path algorithm.
    // Runs until unassigned column added to shortest path tree.
    for (lap_col j = 0; j < dim; ++j) {
      d[j] = input_cost(free_row, j) - v[j];
      predecessor[j] = free_row;
      col_list[j] = j;        // Init column list.
    }
    
    cost min = 0;
    lap_col low = 0; // Columns in 0..low-1 are ready, now none.
    lap_col up = 0;  // Columns in low..up-1 are to be scanned for current minimum, now none.
    // Columns in up..dim-1 are to be considered later to find new minimum;
    // at this stage the list simply contains all columns.
    unassignedfound = false;
    
    do {
      if (up == low) { // No more columns to be scanned for current minimum.
        last = low - 1;
        
        // Scan columns for up..dim-1 to find all indices for which new minimum occurs.
        // Store these indices between low..up-1 (increasing up).
        min = d[col_list[up++]];
        
        for (lap_dim k = up; k < dim; ++k) {
          const lap_col j = col_list[k];
          const cost h = d[j];
          if (h <= min) {
            if (h < min) {   // New minimum.
              up = low;      // Restart list at index low.
              min = h;
            }
            // New index with same minimum, put on undex up, and extend list.
            col_list[k] = col_list[up];
            col_list[up++] = j;
          }
        }
        // Check if any of the minimum columns happens to be unassigned.
        // If so, we have an augmenting path right away.
        for (lap_dim k = low; k < up; ++k) {
          if (colsol[col_list[k]] < 0) {
            endofpath = col_list[k];
            unassignedfound = true;
            break;
          }
        }
      }
      
      if (!unassignedfound) {
        // Update 'distances' between free_row and all unscanned columns,
        // via next scanned column.
        j1 = col_list[low++];
        i = colsol[j1];
        const cost h = input_cost(i, j1) - v[j1] - min;
        
        for (lap_dim k = up; k < dim; ++k) {
          const lap_col j = col_list[k];
          cost v2 = input_cost(i, j) - v[j] - h;
          if (v2 < d[j]) {
            predecessor[j] = i;
            if (v2 == min) { // New column found at same minimum value
              if (colsol[j] < 0) {
                // If unassigned, shortest augmenting path is complete.
                endofpath = j;
                unassignedfound = true;
                break;
              } else {
              // Else add to list to be scanned right away.
                col_list[k] = col_list[up];
                col_list[up++] = j;
              }
            }
            d[j] = v2; // <MS: Unintended>
          }
        }
      }
    } while (!unassignedfound);
    
    // Update column prices.
    for(lap_dim k = 0; k <= last; ++k) {
      j1 = col_list[k];
      v[j1] += d[j1] - min;
    }
    
    // Reset row and column assignments along the alternating path.
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
    const lap_dim j = rowsol[i];
    const cost element_cost = input_cost(i, j);
    u[i] = element_cost - v[j];
    lapcost += element_cost;
  }

  return lapcost;
}
