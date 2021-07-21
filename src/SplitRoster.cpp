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
    splits.emplace_back(xi, n_tips);
    splits[i].quicksort();
  }
  n_bins = splits[0].n_bins;
  
  const int16 max_splits = n_tips - 3;
  roster_tree = std::make_unique<int32[]>(n_trees * max_splits);
  roster_split = std::make_unique<int16[]>(n_trees * max_splits);
  in_split = std::make_unique<int16[]>(n_trees * max_splits);
  out_split = std::make_unique<int16[]>(n_trees * max_splits);
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
  
  
  // Populate leaves of tournament tree
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
  
  roster_len = 0;
  roster_tree[0] = winners[0];
  roster_split[0] = which_split[winners[0]];
  in_split[0] = splits[winners[0]].in_split[roster_split[0]];
  out_split[0] = n_tips - in_split[0];
  roster_hits[0] = 1;
  
  for (; ; ) {
    const int32 winner = winners[0];
    --(which_split[winner]);
    
    int32 i = winner;
    do {
      i = (i - 1) / 2;
      play_game(&i, winners, losers, which_split);
    } while (i);
    
    if (which_split[winners[0]] < 0) break;
    
    push(winners[0], which_split);
  }
  
  ++roster_len; // Point to first position after array
  score.reserve(roster_len * roster_len);
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
    child1 = *node * 2 + 1,
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
                   new_split,
                   splits[roster_tree[roster_len]].state,
                   roster_split[roster_len])) {
    ++roster_hits[roster_len];
  } else {
    ++roster_len;
    roster_tree[roster_len] = tree;
    roster_split[roster_len] = new_split;
    in_split[roster_len] = splits[tree].in_split[new_split];
    out_split[roster_len] = n_tips - in_split[roster_len];
    roster_hits[roster_len] = 1;
  }
  index[tree][new_split] = roster_len;
}

#define SCORE(a, b) score[(a) * roster_len + (b)]

#define SPLIT(i) splits[roster_tree[(i)]].state[roster_split[(i)]]

void SplitRoster::mutual_clustering() {
  const cost max_score = BIG;
  
  for (int16 ai = 0; ai != roster_len - 1; ++ai) {
    const int16
      na = in_split[ai],
      nA = out_split[ai]
    ;
    
    SCORE(ai, ai) = ic_matching(na, nA, n_tips);
    
    for (int16 bi = ai + 1; bi != roster_len; ++bi) {
      
      // x divides tips into a|A; y divides tips into b|B
      int16 a_and_b = 0;
      for (int16 bin = n_bins; bin--; ) {
        a_and_b += count_bits(SPLIT(ai)[bin] & SPLIT(bi)[bin]);
      }
      
      const int16
        nb = in_split[bi],
        nB = out_split[bi],
        a_and_B = na - a_and_b,
        A_and_b = nb - a_and_b,
        A_and_B = nA - A_and_b
      ;
      
      if (a_and_b == A_and_b &&
          a_and_b == a_and_B &&
          a_and_b == A_and_B) {
        SCORE(ai, bi) = max_score; // Don't risk rounding error
      } else {
        SCORE(ai, bi) = max_score -
          // Division by n_tips converts n(A&B) to P(A&B) for each ic_element
          cost(max_score * ((
            // 0 < Sum of IC_elements <= n_tips
            ic_element(a_and_b, na, nb, n_tips) +
            ic_element(a_and_B, na, nB, n_tips) +
            ic_element(A_and_b, nA, nb, n_tips) +
            ic_element(A_and_B, nA, nB, n_tips)
          ) / n_tips)
        );
      }
    }
  }
}

NumericVector SplitRoster::score_pairs() {
  int32 
    i = 0,
    entry = 0,
    n_scores = n_trees * (n_trees + 1) / 2;
  ;
  const cost max_score = BIG;
  //score = std::vector<double> (n_scores, 0.0);
  NumericVector ret(n_scores);
  for (i = 0; ; ++i) {
    const int16 i_splits = splits[i].n_splits;
    // Calculate tree's similarity to self
    for (int16 sp = 0; sp != i_splits; ++sp) {
      const int32 sp_i = index[i][sp];
      ret[entry] += SCORE(sp_i, sp_i);
    }
    
    if (i == n_trees) break;
    entry++;
    for (int32 j = i + 1; j != n_trees; ++j) {
      const int16
        j_splits = splits[j].n_splits,
        most_splits = i_splits > j_splits ? i_splits : j_splits
      ;
      
      
      // TODO declare lap_score once at maximum size n_tips - 3 and retain.
      cost** lap_score = new cost*[most_splits];
      for (int16 i = most_splits; i--; ) lap_score[i] = new cost[most_splits];
      lap_col *rowsol = new lap_col[most_splits];
      lap_row *colsol = new lap_row[most_splits];
      cost *u = new cost[most_splits], *v = new cost[most_splits];
      
      
      
      for (int16 ai = i_splits; ai--; ) {
        for (int16 bi = j_splits; bi--; ) {
          lap_score[ai][bi] = SCORE(ai, bi);
        }
      }
      
      for (int16 ai = i_splits; ai < most_splits; ++ai) {
        for (int16 bi = 0; bi != most_splits; ++bi) {
          lap_score[ai][bi] = max_score;
        }
      }
      
      ret[entry] = double(
        (max_score * most_splits) -
          lap(most_splits, lap_score, rowsol, colsol, u, v)
      ) / max_score;
      
      for (int16 i = most_splits; i--; ) delete[] lap_score[i];
      delete[] colsol; delete[] u; delete[] v; delete[] lap_score;
      delete[] rowsol;
      
      entry++;
    }
  }
  return ret;
}
