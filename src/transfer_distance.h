#ifndef TREEDIST_TRANSFER_DISTANCE_H
#define TREEDIST_TRANSFER_DISTANCE_H

#include <TreeTools/SplitList.h>

// Transfer dissimilarity between two trees; thread-safe.
double transfer_dissimilarity_score(
    const TreeTools::SplitList& sl_x, const TreeTools::SplitList& sl_y,
    int n_tip, int n_bins, bool scale
);

#endif
