#ifndef TREEDIST_SPLIT_COUNTS_H_
#define TREEDIST_SPLIT_COUNTS_H_

#include <cstdint>

// Splits of a pair of trees, each obtained from a parent tree by removing
// one leaf, described through the parent's splits.  The trees have `n_split`
// splits each.  Child split i of the first tree restricts parent split
// row[i]; it has row_in[i] leaves, and row_has[i] is 1 if the parent split
// contained the removed leaf, else 0.  col, col_in and col_has describe the
// second tree likewise.  overlap[r * stride + c] counts the leaves shared by
// parent splits r and c, so the corresponding child splits share
//   overlap[r * stride + c] - (row_has & col_has).
struct SplitCounts {
  int n_tip; // Leaves in each child tree
  int n_split;
  const int16_t *overlap;
  int stride;
  const int *row;
  const int *row_in;
  const int16_t *row_has;
  const int *col;
  const int *col_in;
  const int16_t *col_has;
};

#endif
