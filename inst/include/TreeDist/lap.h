#ifndef TREEDIST_LAP_H_
#define TREEDIST_LAP_H_

// LAP (Linear Assignment Problem) — Jonker–Volgenant declarations.
//
// Provides the lap() function signature.  The implementation lives in
// lap_impl.h, guarded by TREEDIST_LAP_IMPLEMENTATION.  Downstream
// packages should include lap_impl.h in exactly one translation unit
// with the guard defined.

#include "lap_scratch.h"

namespace TreeDist {

  // Primary overload: caller supplies pre-allocated scratch.
  extern cost lap(lap_row dim,
                  CostMatrix& input_cost,
                  std::vector<lap_col>& rowsol,
                  std::vector<lap_row>& colsol,
                  bool allow_interrupt,
                  LapScratch& scratch);

  // Convenience overload: creates a temporary scratch.
  inline cost lap(lap_row dim,
                  CostMatrix& input_cost,
                  std::vector<lap_col>& rowsol,
                  std::vector<lap_row>& colsol,
                  bool allow_interrupt = true) {
    LapScratch scratch;
    return lap(dim, input_cost, rowsol, colsol, allow_interrupt, scratch);
  }

} // namespace TreeDist

#endif // TREEDIST_LAP_H_
