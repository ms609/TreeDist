#ifndef TREEDIST_TYPES_H_
#define TREEDIST_TYPES_H_

// TreeDist public type definitions.
//
// Rcpp-free: uses only standard headers plus TreeTools (for SL_MAX_SPLITS).
// Downstream packages include via:  #include <TreeDist/types.h>

#include <cstdint>
#include <limits>
#include <TreeTools/SplitList.h>   // SL_MAX_SPLITS, splitbit

namespace TreeDist {

  using int16 = int_fast16_t;
  using int32 = int_fast32_t;
  using cost  = int_fast64_t;

  // Canonical type for split/tip/bin counters.
  using split_int = int32;

  using lap_dim = int;
  using lap_row = lap_dim;
  using lap_col = lap_dim;

  constexpr cost BIG =
    (std::numeric_limits<cost>::max)() / SL_MAX_SPLITS;

  constexpr cost ROUND_PRECISION = 2048 * 2048;

  constexpr std::size_t BLOCK_SIZE = 16;

} // namespace TreeDist

#endif // TREEDIST_TYPES_H_
