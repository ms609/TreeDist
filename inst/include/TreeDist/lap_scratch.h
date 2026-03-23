#ifndef TREEDIST_LAP_SCRATCH_H_
#define TREEDIST_LAP_SCRATCH_H_

// Reusable heap storage for one thread's LAPJV workspace.
//
// Pass to the scratch-taking overload of lap() to eliminate per-call
// std::vector allocations.  In a serial context construct once before the
// loop; in an OpenMP context allocate one per thread and index by
// omp_get_thread_num().  All vectors are grown lazily; never shrunk.

#include "cost_matrix.h"

namespace TreeDist {

struct LapScratch {
  std::vector<cost>    v;
  std::vector<lap_col> matches;
  std::vector<lap_row> freeunassigned;
  std::vector<lap_col> col_list;
  std::vector<cost>    d;
  std::vector<lap_row> predecessor;
  std::vector<lap_col> rowsol;
  std::vector<lap_row> colsol;
  CostMatrix score_pool;
  CostMatrix small_pool;

  void ensure(int dim) noexcept {
    const int padded = static_cast<int>(
      ((static_cast<std::size_t>(dim) + BLOCK_SIZE - 1)
       / BLOCK_SIZE) * BLOCK_SIZE);
    if (static_cast<int>(v.size())              < padded) v.resize(padded);
    if (static_cast<int>(matches.size())        < dim)    matches.resize(dim);
    if (static_cast<int>(freeunassigned.size()) < dim)    freeunassigned.resize(dim);
    if (static_cast<int>(col_list.size())       < dim)    col_list.resize(dim);
    if (static_cast<int>(d.size())              < dim)    d.resize(dim);
    if (static_cast<int>(predecessor.size())    < dim)    predecessor.resize(dim);
    if (static_cast<int>(rowsol.size())         < dim)    rowsol.resize(dim);
    if (static_cast<int>(colsol.size())         < dim)    colsol.resize(dim);
  }
};

} // namespace TreeDist

#endif // TREEDIST_LAP_SCRATCH_H_
