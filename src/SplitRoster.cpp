#include <Rcpp.h>
#include "ints.hpp"
#include "tree_distances.hpp"
#include "SplitList.hpp"
#include "SplitRoster.hpp"
using namespace Rcpp;

inline bool SplitRoster::game_result(
    const splitbit (&a)[MAX_SPLITS][MAX_BINS], const int16 split_a,
    const splitbit (&b)[MAX_SPLITS][MAX_BINS], const int16 split_b) {
  for (int16 bin = 0; bin != n_bins; ++bin) {
    if (a[split_a][bin] > b[split_b][bin]) {
      return true;
    }
  }
  return false;
}

inline bool SplitRoster::splits_equal(
    const splitbit (&a)[MAX_SPLITS][MAX_BINS], const int16 split_a,
    const splitbit (&b)[MAX_SPLITS][MAX_BINS], const int16 split_b) {
  for (int16 bin = n_bins; bin--; ) {
    if (a[split_a][bin] != b[split_b][bin]) {
      return false;
    }
  }
  return true;
}


inline int32 SplitRoster::game_winner(
    const int32 *node,
    std::unique_ptr<int32[]> &which_tree,
    std::unique_ptr<int16[]> &which_split) {
  const int32
    child1 = *node * 2,
    child2 = child1 + 1
  ;
  
  bool child1_greater;
  if (which_split[child1] < 0) {
    child1_greater = false;
  } else if (which_split[child2] < 0) {
    child1_greater = true;
  } else {
    child1_greater = game_result(
      splits[which_tree[child1]].state, which_split[child1],
      splits[which_tree[child2]].state, which_split[child2]));
  }
  
  const int32 loser = child1_greater ? child2 : child1;
  which_tree[*node] = which_tree[loser];
  which_split[*node] = which_split[loser];
  return child1_greater ? child1 : child2;
}

inline void SplitRoster::push(
    const int32 node,
    std::unique_ptr<int32[]> &which_tree,
    std::unique_ptr<int16[]> &which_split) {
  const int32 new_tree = which_tree[node];
  const int16 new_split = which_split[node];
  if (splits_equal(splits[which_tree[new_tree]].state,
                   which_split[new_split],
                   splits[roster_tree[roster_pos]].state,
                   roster_split[roster_pos])) {
    ++roster_hits[roster_pos];
  } else {
    ++roster_pos;
    roster_tree[roster_pos] = new_tree;
    roster_split[roster_pos] = new_split;
    roster_size[roster_pos] = splits[which_tree[new_tree]].in_split;
    roster_hits[roster_pos] = 1;
  }
  index[new_tree][new_split] = roster_pos;
}

// [[Rcpp::export]]
SplitRoster::SplitRoster(const List x, const IntegerVector nTip) {
  n_tips = nTip[0];
  n_trees = x.length();
  splits = std::make_unique<SplitList[]>(n_trees);
  
  for (int32 i = n_trees; i--; ) {
    const RawMatrix a = x[i];
    splits[i] = SplitList(a);
    splits[i].quicksort();
  }
  n_bins = splits[0].n_bins;
  
  const int16 max_splits = n_tips - 3;
  roster_tree = std::make_unique<int32[]>(n_trees * max_splits);
  roster_split = std::make_unique<int16[]>(n_trees * max_splits);
  roster_size = std::make_unique<int16[]>(n_trees * max_splits);
  roster_hits = std::make_unique<int32[]>(n_trees * max_splits);
  index = std::make_unique<std::unique_ptr<int16>[]>(n_trees);
  for (int32 i = n_trees; i--; ) {
    index[i] = std::make_unique<int16>(max_splits);
  }
  
  // Populate roster using k-way merge with tournament tree
  const int32
    tournament_games = n_trees - 1,
    tournament_nodes = n_trees + tournament_games
  ;
  
  auto which_tree = std::make_unique<int16[]>(tournament_nodes);
  auto which_split = std::make_unique<int16[]>(tournament_nodes);
  auto winner_index = std::make_unique<int16[]>(tournament_nodes);
  
  
  // Populate children of tree
  for (int32 i = 0; i != n_trees; ++i) {
    const int32 this_node = i + tournament_games;
    which_tree[this_node] = i;
    which_split[this_node] = splits[i].n_splits - 1;
    winner_index[this_node] = this_node;
  }
  for (int32 i = n_trees + tournament_games; i != tournament_nodes; ++i) {
    which_split[i] = -1;
  }
  
  // Initial games
  for (int32 i = tournament_games; i--; ) {
    // TODO duplicate function as void to avoid overhead of unused return
    winner_index[i] = game_winner(&i, which_tree, which_split);
  }
  
  roster_pos = 0;
  roster_tree[0] = which_tree[winner_index[0]];
  roster_split[0] = which_split[winner_index[0]];
  roster_hits[0] = 1;
  
  for (; ; ) {
    const int32 winner = winner_index[0];
    
    
    int32 i = winner;
    do {
      i /= 2;
      winner_index[i] = game_winner(&i, which_tree, which_split);
    } while (i);
    
    push(winner_index[0], which_tree, which_split);
  }
  
  
  
  
}
