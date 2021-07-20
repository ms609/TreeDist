#include <Rcpp.h>
#include "ints.hpp"
#include "tree_distances.hpp"
#include "SplitList.hpp"
#include "SplitRoster.hpp"
using namespace Rcpp;

SplitRoster::SplitRoster(const List x, const IntegerVector nTip) {
  n_tips = nTip[0];
  n_trees = x.length();
  splits.reserve(n_trees);
  
  for (int32 i = 0; i != n_trees; ++i) {
    const RawMatrix xi = x[i];
    splits.emplace_back(xi);
    splits[i].quicksort();
  }
  n_bins = splits[0].n_bins;
  
  const int16 max_splits = n_tips - 3;
  roster_tree = std::make_unique<int32[]>(n_trees * max_splits);
  roster_split = std::make_unique<int16[]>(n_trees * max_splits);
  roster_size = std::make_unique<int16[]>(n_trees * max_splits);
  roster_hits = std::make_unique<int32[]>(n_trees * max_splits);
  index.reserve(n_trees);
  
  // Populate roster using k-way merge with tournament tree
  const int32
    tournament_games = n_trees - 1,
      tournament_nodes = n_trees + tournament_games
    ;
  
  auto winners = std::make_unique<int32[]>(tournament_nodes);
  auto losers = std::make_unique<int32[]>(tournament_nodes);
  auto which_split = std::make_unique<int16[]>(n_trees);
  
  
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
    play_game(&i, winners, losers, which_split);
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
      play_game(&i, winners, losers, which_split);
    } while (i);
    
    if (which_split[winners[0]] < 0) break;
    
    push(winners[0], which_split);
  }
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
      splits[tree2].state, which_split[tree2]);
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
    const int32 tree,
    std::unique_ptr<int16[]> &which_split) {
  const int16 new_split = which_split[tree];
  if (splits_equal(splits[tree].state,
                   which_split[new_split],
                   splits[roster_tree[roster_pos]].state,
                   roster_split[roster_pos])) {
    ++roster_hits[roster_pos];
  } else {
    ++roster_pos;
    roster_tree[roster_pos] = tree;
    roster_split[roster_pos] = new_split;
    roster_size[roster_pos] = splits[tree].in_split[new_split];
    roster_hits[roster_pos] = 1;
  }
  index[tree][new_split] = roster_pos;
}

  
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
    play_game(&i, winners, losers, which_split);
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
      play_game(&i, winners, losers, which_split);
    } while (i);
    
    if (which_split[winners[0]] < 0) break;
    
    push(winners[0], which_split);
  }
}
