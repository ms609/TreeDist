#ifndef _TREEDIST_LAP_H
#define _TREEDIST_LAP_H

// LAP (Linear Assignment Problem) infrastructure.
//
// This header is intentionally kept separate from the tree-distance scoring
// functions in tree_distances.h.  Changes to scoring code (add_ic_element,
// one_overlap, spi_overlap, etc.) must not affect the compilation of lap.cpp,
// because the LAP hot loops are sensitive to instruction alignment — even
// unrelated inline functions in the same translation unit can shift code
// layout enough to cause measurable regressions.
//
// Public types, CostMatrix, LapScratch, and lap() declarations live in
// inst/include/TreeDist/ headers (shared with downstream LinkingTo users).
// This wrapper re-exports them to global scope for backward compatibility.

#include <TreeDist/lap.h>     // types, CostMatrix, LapScratch, lap()

#include "ints.h"             // grf_match, uint32 (TreeDist-internal types)

// Re-export shared types and classes to global scope.
using TreeDist::cost;
using TreeDist::lap_dim;
using TreeDist::lap_row;
using TreeDist::lap_col;
using TreeDist::BIG;
using TreeDist::ROUND_PRECISION;
using TreeDist::BLOCK_SIZE;
using TreeDist::CostMatrix;
using TreeDist::LapScratch;
using TreeDist::lap;

using cost_matrix = TreeDist::CostMatrix;

#endif
