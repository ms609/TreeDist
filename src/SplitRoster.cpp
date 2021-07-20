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


inline void SplitRoster::play_game(
    const int32 *node,
    std::unique_ptr<int32[]> &winners,
    std::unique_ptr<int32[]> &losers,
    std::unique_ptr<int16[]> &which_split) {
  const int32
    child1 = *node * 2,
    child2 = child1 + 1,
    tree1 = winners[child1],
    tree2 = winners[child2]
  ;
  
  bool child1_greater;
  if (tree1 < 0 || which_split[tree1] < 0) {
    child1_greater = false;
  } else if (tree2 < 0 || which_split[tree2] < 0) {
    child1_greater = true;
  } else {
    child1_greater = game_result(
      splits[tree1].state, which_split[tree1],
      splits[tree2].state, which_split[tree2]));
  }
  
  if (child1_greater) {
    winners[*node] = tree1;
    losers[*node] = tree2;
  } else {
    winners[*node] = tree2;
    losers[*node] = tree1;
  }
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
  
  auto which_split = std::make_unique<int16[]>(n_trees);
  auto winners = std::make_unique<int16[]>(tournament_nodes);
  auto losers = std::make_unique<int16[]>(tournament_nodes);
  
  
  // Populate children of tree
  for (int32 i = 0; i != n_trees; ++i) {
    which_split[i] = splits[i].n_splits - 1;
    winners[i + tournament_games] = i;
  }
  for (int32 i = n_trees + tournament_games; i != tournament_nodes; ++i) {
    which_split[i] = -1;
    winners[i] = -1;
  }
  
  // Initial games
  for (int32 i = tournament_games; i--; ) {
    play_game(&i, which_split, winners, losers);
  }
  
  roster_pos = 0;
  roster_tree[0] = winners[0];
  roster_split[0] = which_split[winners[0]];
  roster_hits[0] = 1;
  
  for (; ; ) {
    const int32 winner = winners[0];
    --which_split[winner];
    
    int32 i = winner;
    do {
      i /= 2;
      play_game(&i, which_tree, which_split, winners, losers);
    } while (i);
    
    if (which_split[winners[0]] < 0) break;
    
    push(winners[0], which_tree, which_split);
  }
  
  
  
  
}
