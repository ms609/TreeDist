#ifndef _TREEDIST_SPLITROSTER_HPP
#define _TREEDIST_SPLITROSTER_HPP

#include <stdint.h>
#include <Rcpp.h>
#include "ints.hpp"
#include "tree_distances.hpp"
#include "SplitList.hpp"

using namespace Rcpp;

class SplitRoster {
  bool game_result(const splitbit (&a)[MAX_SPLITS][MAX_BINS], const int16 split_a,
                   const splitbit (&b)[MAX_SPLITS][MAX_BINS], const int16 split_b);
  int32 game_winner(
      const int32 *node,
      std::unique_ptr<int32[]> &which_tree,
      std::unique_ptr<int16[]> &which_split);
  void push();
  std::unique_ptr<SplitList[]> splits;
  std::unique_ptr<int32[]> roster_tree;
  std::unique_ptr<int16[]> roster_split;
  std::unique_ptr<int16[]> roster_size;
  std::unique_ptr<int32[]> roster_hits;
  std::unique_ptr<std::unique_ptr<int16>[]> index;
  int16 n_bins;
  int32 roster_pos;
public:
  int16 n_tips, n_splits;
  int32 n_trees;
};

#endif
