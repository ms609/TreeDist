#ifndef _TREEDIST_SPLITROSTER_HPP
#define _TREEDIST_SPLITROSTER_HPP

#include <vector>
#include <stdint.h>
#include <Rcpp.h>
#include "ints.hpp"
#include "tree_distances.hpp"
#include "SplitList.hpp"

using namespace Rcpp;

class SplitRoster {
  int16 n_bins;
  int32 roster_pos;
  
  std::unique_ptr<SplitList[]> splits;
  std::unique_ptr<int32[]> roster_tree;
  std::unique_ptr<int16[]> roster_split;
  std::unique_ptr<int16[]> roster_size;
  std::unique_ptr<int32[]> roster_hits;
  std::vector<std::array<int32, MAX_SPLITS>[]> index;
  
  bool splits_equal(
      const splitbit (&a)[MAX_SPLITS][MAX_BINS], const int16 split_a,
      const splitbit (&b)[MAX_SPLITS][MAX_BINS], const int16 split_b);
  bool game_result(
      const splitbit (&a)[MAX_SPLITS][MAX_BINS], const int16 split_a,
      const splitbit (&b)[MAX_SPLITS][MAX_BINS], const int16 split_b);
  void play_game(
      const int32 *node,
      std::unique_ptr<int32[]> &winners,
      std::unique_ptr<int32[]> &losers,
      std::unique_ptr<int16[]> &which_split);
  void push(
      const int32 tree,
      std::unique_ptr<int16[]> &which_split);
  
public:
  int16 n_tips, n_splits;
  int32 n_trees;
  
  SplitRoster(const List x, const IntegerVector nTip);
};

#endif
