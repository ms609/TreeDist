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
  in_split = std::make_unique<int16[]>(n_trees * max_splits);
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
  
  roster_len = 0;
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
                   splits[roster_tree[roster_len]].state,
                   roster_split[roster_len])) {
    ++roster_hits[roster_len];
  } else {
    ++roster_len;
    roster_tree[roster_len] = tree;
    roster_split[roster_len] = new_split;
    in_split[roster_len] = splits[tree].in_split[new_split];
    roster_hits[roster_len] = 1;
  }
  index[tree][new_split] = roster_len;
}

#define SCORE(a, b) score[(a) * roster_len + (b)]

#define SPLIT(i) splits[roster_tree[(i)]].state[roster_split[(i)]]

void SplitRoster::mutual_clustering() {
  const cost max_score = BIG;
  double exact_match_score = 0;
  int16 exact_matches = 0;
  
  for (int16 ai = 0; ai != roster_len - 1; ++ai) {
    const int16
      na = in_split[ai],
      nA = n_tips - na
    ;
    
    SCORE(ai, ai) = ic_matching(na, nA, n_tips);
    
    for (int16 bi = ai + 1; bi != roster_len; ++bi) {
      
      // x divides tips into a|A; y divides tips into b|B
      int16 a_and_b = 0;
      for (int16 bin = n_bins; --bin; ) {
        a_and_b += count_bits(SPLIT(ai)[bin] & SPLIT(bi)[bin]);
      }
      
      const int16
        nb = in_split[bi],
        nB = n_tips - nb,
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


// [[Rcpp::export]]
NumericVector SplitRoster::score_pairs() {
  int32 i = 0, entry = 0;
  for (i = 0; ; ++i) {
    // Calculate tree's similarity to self
    
    
    if (i == n_trees) break;
    entry++;
    for (int32 j = i + 1; j != n_trees; ++j) {
      
      entry++;
    }
  }
  const SplitList a(x), b(y);
  const bool a_has_more_splits = (a.n_splits > b.n_splits);
  const int16
    most_splits = a_has_more_splits ? a.n_splits : b.n_splits,
      a_extra_splits = a_has_more_splits ? most_splits - b.n_splits : 0,
      b_extra_splits = a_has_more_splits ? 0 : most_splits - a.n_splits,
      n_tips = nTip[0]
  ;
  const cost max_score = BIG;
  
  cost** score = new cost*[most_splits];
  for (int16 i = most_splits; i--; ) score[i] = new cost[most_splits];
  double exact_match_score = 0;
  int16 exact_matches = 0;
  // NumericVector zero-initializes [so does make_unique]
  // match will have one added to it so numbering follows R; hence 0 = UNMATCHED
  NumericVector a_match(a.n_splits);
  std::unique_ptr<int16[]> b_match = std::make_unique<int16[]>(b.n_splits);
  
  for (int16 ai = 0; ai != a.n_splits; ++ai) {
    if (a_match[ai]) continue;
    const int16
      na = a.in_split[ai],
                     nA = n_tips - na
      ;
    
    for (int16 bi = 0; bi != b.n_splits; ++bi) {
      
      // x divides tips into a|A; y divides tips into b|B
      int16 a_and_b = 0;
      for (int16 bin = 0; bin != a.n_bins; ++bin) {
        a_and_b += count_bits(a.state[ai][bin] & b.state[bi][bin]);
      }
      
      const int16
        nb = b.in_split[bi],
                       nB = n_tips - nb,
                       a_and_B = na - a_and_b,
                       A_and_b = nb - a_and_b,
                       A_and_B = nA - A_and_b
        ;
      
      if ((!a_and_B && !A_and_b) ||
          (!a_and_b && !A_and_B)) {
        exact_match_score += ic_matching(na, nA, n_tips);
        exact_matches++;
        a_match[ai] = bi + 1;
        b_match[bi] = ai + 1;
        break;
      } else if (a_and_b == A_and_b &&
        a_and_b == a_and_B &&
        a_and_b == A_and_B) {
        score[ai][bi] = max_score; // Don't risk rounding error
      } else {
        score[ai][bi] = max_score -
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
    for (int16 bi = b.n_splits; bi < most_splits; ++bi) {
      score[ai][bi] = max_score;
    }
  }
  if (exact_matches == b.n_splits || exact_matches == a.n_splits) {
    for (int16 i = most_splits; i--; ) delete[] score[i];
    delete[] score;
    
    return List::create(
      Named("score") = NumericVector::create(exact_match_score / n_tips),
      _["matching"] = a_match);
  }
  
  
  const int16 lap_dim = most_splits - exact_matches;
  lap_col *rowsol = new lap_col[lap_dim];
  lap_row *colsol = new lap_row[lap_dim];
  cost *u = new cost[lap_dim], *v = new cost[lap_dim];
  
  if (exact_matches) {
    int16 a_pos = 0;
    for (int16 ai = 0; ai != a.n_splits; ++ai) {
      if (a_match[ai]) continue;
      int16 b_pos = 0;
      for (int16 bi = 0; bi != b.n_splits; ++bi) {
        if (b_match[bi]) continue;
        score[a_pos][b_pos] = score[ai][bi];
        b_pos++;
      }
      for (int16 bi = lap_dim - a_extra_splits; bi < lap_dim; ++bi) {
        score[a_pos][bi] = max_score;
      }
      a_pos++;
    }
    for (int16 ai = lap_dim - b_extra_splits; ai < lap_dim; ++ai) {
      for (int16 bi = 0; bi != lap_dim; ++bi) {
        score[ai][bi] = max_score;
      }
    }
    
    const double lap_score = 
      double((max_score * lap_dim) - lap(lap_dim, score, rowsol, colsol, u, v))
      / max_score;
    NumericVector final_score = 
    NumericVector::create(lap_score + (exact_match_score / n_tips));
    
    for (int16 i = most_splits; i--; ) delete[] score[i];
    delete[] colsol; delete[] u; delete[] v; delete[] score;
    
    std::unique_ptr<int16[]> lap_decode = std::make_unique<int16[]>(lap_dim);
    int16 fuzzy_match = 0;
    for (int16 bi = 0; bi != b.n_splits; ++bi) {
      if (!b_match[bi]) {
        lap_decode[fuzzy_match++] = bi + 1;
      } else {
      }
    }
    
    fuzzy_match = 0;
    IntegerVector final_matching(a.n_splits);
    for (int16 i = 0; i != a.n_splits; i++) {
      if (a_match[i]) {
        // Rcout << "a" << (1+i) << " exactly matches b" << a_match[i]<< "\n";
        final_matching[i] = a_match[i];
      } else {
        const int16 this_sol = rowsol[fuzzy_match++];
        // Rcout << "a"<<(1+i) << " fuzzily matches rowsol[" << this_sol <<"] == "
        //       << rowsol[this_sol] << "; ";
        if (rowsol[this_sol] >= lap_dim - a_extra_splits) {
          // Rcout << " unmatched (NA)\n";
          final_matching[i] = NA_INTEGER;
        } else {
          // Rcout << " matched with b" << lap_decode[rowsol[this_sol]] <<".\n";
          final_matching[i] = lap_decode[rowsol[this_sol]];
        }
      }
      // Rcout << " ";
      // if (final_matching[i] > 0) Rcout << final_matching[i]; else Rcout << "NA";
    }
    
    delete[] rowsol;
    
    return List::create(Named("score") = final_score,
                        _["matching"] = final_matching);
  } else {
    for (int16 ai = a.n_splits; ai < most_splits; ++ai) {
      for (int16 bi = 0; bi != most_splits; ++bi) {
        score[ai][bi] = max_score;
      }
    }
    
    const double lap_score = double(
      (max_score * lap_dim) -
        lap(most_splits, score, rowsol, colsol, u, v)
    ) / max_score;
    NumericVector final_score = NumericVector::create(lap_score);
    
    for (int16 i = most_splits; i--; ) delete[] score[i];
    delete[] colsol; delete[] u; delete[] v; delete[] score;
    
    IntegerVector final_matching (a.n_splits);
    for (int16 i = a.n_splits; i--; ) {
      final_matching[i] = (rowsol[i] < b.n_splits) ? rowsol[i] + 1 : NA_INTEGER;
    }
    
    delete[] rowsol;
    
    return List::create(Named("score") = final_score,
                        _["matching"] = final_matching);
  }
}
