#include <cmath>
#include <Rcpp.h>
#include "tree_distances.h"
#include "SplitList.h"

using namespace Rcpp;

// [[Rcpp::export]]
List cpp_robinson_foulds_distance (RawMatrix x, RawMatrix y, 
                                   IntegerVector nTip) {
  if (x.cols() != y.cols()) {
    throw std::invalid_argument("Input splits must address same number of tips.");
  }
  SplitList a(x), b(y);
  const int16 most_splits = (a.n_splits > b.n_splits) ? a.n_splits : b.n_splits,
              last_bin = a.n_bins - 1,
              n_tips = nTip[0],
              unset_tips = (n_tips % BIN_SIZE) ? BIN_SIZE - n_tips % BIN_SIZE : 0;
  const splitbit unset_mask = ALL_ONES >> unset_tips;
  cost score = 0;
  
  rf_match matching (most_splits);
  for (int16 i = 0; i != most_splits; i++) matching[i] = NA_INTEGER;
  
  splitbit b_complement[MAX_SPLITS][MAX_BINS];
  for (int16 i = 0; i != b.n_splits; i++) {
    for (int16 bin = 0; bin != last_bin; bin++) {
      b_complement[i][bin] = ~b.state[i][bin];
    }
    b_complement[i][last_bin] = b.state[i][last_bin] ^ unset_mask;
  }
  
  for (int16 ai = 0; ai != a.n_splits; ai++) {
    for (int16 bi = 0; bi != b.n_splits; bi++) {
    
      bool all_match = true, all_complement = true;
    
      for (int16 bin = 0; bin != a.n_bins; bin++) {
        if ((a.state[ai][bin] != b.state[bi][bin])) {
          all_match = false;
          break;
        }
      }
      if (!all_match) {
        for (int16 bin = 0; bin != a.n_bins; bin++) {
          if ((a.state[ai][bin] != b_complement[bi][bin])) {
            all_complement = false;
            break;
          }
        }
      }
      if (all_match || all_complement) {
        ++score;
        matching[ai] = bi + 1;
        break; /* Only one match possible per split */
      }
    }
  }
  score = cost (a.n_splits + b.n_splits) - score - score;
  
  List ret = List::create(Named("score") = Rcpp::wrap(score),
                          _["matching"] = Rcpp::wrap(matching));
  
  return (ret);
}

rf_match cpp_robinson_foulds_matching (
    raw_vector x, raw_vector y,
    const int16 n_cols,
    const int16 n_tips) {
  
  SplitList a(x, n_cols), b(y, n_cols);
  const int16 most_splits = (a.n_splits > b.n_splits) ? a.n_splits : b.n_splits,
    last_bin = a.n_bins - 1,
    unset_tips = (n_tips % BIN_SIZE) ? BIN_SIZE - n_tips % BIN_SIZE : 0;
  const splitbit unset_mask = ALL_ONES >> unset_tips;
  cost score = 0;
  
  rf_match matching (most_splits);
  for (int16 i = 0; i != most_splits; i++) matching[i] = NA_INTEGER;
  
  splitbit b_complement[MAX_SPLITS][MAX_BINS];
  for (int16 i = 0; i != b.n_splits; i++) {
    for (int16 bin = 0; bin != last_bin; bin++) {
      b_complement[i][bin] = ~b.state[i][bin];
    }
    b_complement[i][last_bin] = b.state[i][last_bin] ^ unset_mask;
  }
  
  for (int16 ai = 0; ai != a.n_splits; ai++) {
    for (int16 bi = 0; bi != b.n_splits; bi++) {
    
      bool all_match = true, all_complement = true;
    
      for (int16 bin = 0; bin != a.n_bins; bin++) {
        if ((a.state[ai][bin] != b.state[bi][bin])) {
          all_match = false;
          break;
        }
      }
      if (!all_match) {
        for (int16 bin = 0; bin != a.n_bins; bin++) {
          if ((a.state[ai][bin] != b_complement[bi][bin])) {
            all_complement = false;
            break;
          }
        }
      }
      if (all_match || all_complement) {
        ++score;
        matching[ai] = bi + 1;
        break; /* Only one match possible per split */
      }
    }
  }
  
  return (matching);
}

// [[Rcpp::export]]
List cpp_robinson_foulds_info (RawMatrix x, RawMatrix y, 
                               IntegerVector nTip) {
  if (x.cols() != y.cols()) {
    throw std::invalid_argument("Input splits must address same number of tips.");
  }
  SplitList a(x), b(y);
  const int16 most_splits = (a.n_splits > b.n_splits) ? a.n_splits : b.n_splits,
              last_bin = a.n_bins - 1,
              n_tips = nTip[0],
              unset_tips = (n_tips % BIN_SIZE) ? BIN_SIZE - n_tips % BIN_SIZE : 0;
  const splitbit unset_mask = ALL_ONES >> unset_tips;
  const double lg2_unrooted_n = lg2_unrooted[n_tips];
  double score = 0;
  
  IntegerVector matching (most_splits);
  for (int16 i = 0; i != most_splits; i++) matching[i] = NA_REAL;
  
  /* Dynamic allocation 20% faster for 105 tips, but VLA not permitted in C11 */
  splitbit b_complement[MAX_SPLITS][MAX_BINS]; 
  for (int16 i = 0; i != b.n_splits; i++) {
    for (int16 bin = 0; bin != last_bin; bin++) {
      b_complement[i][bin] = ~b.state[i][bin];
    }
    b_complement[i][last_bin] = b.state[i][last_bin] ^ unset_mask;
  }
  
  for (int16 ai = 0; ai != a.n_splits; ai++) {
    for (int16 bi = 0; bi != b.n_splits; bi++) {
      
      bool all_match = true, all_complement = true;
      
      for (int16 bin = 0; bin != a.n_bins; bin++) {
        if ((a.state[ai][bin] != b.state[bi][bin])) {
          all_match = false;
          break;
        }
      }
      if (!all_match) {
        for (int16 bin = 0; bin != a.n_bins; bin++) {
          if ((a.state[ai][bin] != b_complement[bi][bin])) {
            all_complement = false;
            break;
          }
        }
      }
      if (all_match || all_complement) {
        int16 leaves_in_split = 0;
        for (int16 bin = 0; bin != a.n_bins; bin++) {
          leaves_in_split += count_bits(a.state[ai][bin]);
        }
        score += lg2_unrooted_n - 
          lg2_trees_matching_split(leaves_in_split, 
                                   n_tips - leaves_in_split);
        matching[ai] = bi + 1;
        break; /* Only one match possible per split */
      }
    }
  }
  
  NumericVector final_score = NumericVector::create(score);
  List ret = List::create(Named("score") = final_score,
                          _["matching"] = matching);
  
  return (ret);
}

// [[Rcpp::export]]
List cpp_matching_split_distance (RawMatrix x, RawMatrix y, 
                                  IntegerVector nTip) {
  if (x.cols() != y.cols()) {
    throw std::invalid_argument("Input splits must address same number of tips.");
  }
  SplitList a(x), b(y);
  const int16 most_splits = (a.n_splits > b.n_splits) ? a.n_splits : b.n_splits,
    split_diff = most_splits - 
      ((a.n_splits > b.n_splits) ? b.n_splits : a.n_splits),
    n_tips = nTip[0],
    half_tips = n_tips / 2;
  const cost max_score = BIG / most_splits;
  
  cost** score = new cost*[most_splits];
  for (int16 i = 0; i < most_splits; i++) score[i] = new cost[most_splits];
  
  for (int16 ai = 0; ai != a.n_splits; ai++) {
    for (int16 bi = 0; bi != b.n_splits; bi++) {
      score[ai][bi] = 0;
      for (int16 bin = 0; bin != a.n_bins; bin++) {
        score[ai][bi] += count_bits(a.state[ai][bin] ^ 
          b.state[bi][bin]);
      }
      if (score[ai][bi] > half_tips) score[ai][bi] = n_tips - score[ai][bi];
    }
    for (int16 bi = b.n_splits; bi < most_splits; bi++) {
      score[ai][bi] = max_score;
    }
  }
  for (int16 ai = a.n_splits; ai < most_splits; ai++) {
    for (int16 bi = 0; bi < most_splits; bi++) {
      score[ai][bi] = max_score;
    }
  }
  
  lap_col *rowsol = new lap_col[most_splits];
  lap_row *colsol = new lap_row[most_splits];
  cost *u = new cost[most_splits], 
                    *v = new cost[most_splits];
  
  NumericVector final_score = NumericVector::create(
    lap(most_splits, score, rowsol, colsol, u, v) - (max_score * split_diff));
  for (int16 i = 0; i < most_splits; i++) delete[] score[i];
  delete[] u; delete[] v; delete[] colsol; delete[] score;
  NumericVector final_matching (most_splits);
  
  for (int16 i = 0; i < most_splits; i++) {
    final_matching[i] = rowsol[i] + 1;
  }
  delete[] rowsol;
  
  List ret = List::create(Named("score") = final_score,
                          _["matching"] = final_matching);
  
  return (ret);
}

// [[Rcpp::export]]
List cpp_jaccard_similarity (RawMatrix x, RawMatrix y,
                             IntegerVector nTip, NumericVector k,
                             LogicalVector arboreal) {
  if (x.cols() != y.cols()) {
    throw std::invalid_argument("Input splits must address same number of tips.");
  }
  SplitList a(x), b(y);
  const int16 most_splits = (a.n_splits > b.n_splits) ? a.n_splits : b.n_splits,
    last_bin = a.n_bins - 1,
    n_tips = nTip[0],
    unset_tips = (n_tips % BIN_SIZE) ? BIN_SIZE - n_tips % BIN_SIZE : 0;
  const cost max_score = BIG;
  const splitbit unset_mask = ALL_ONES >> unset_tips;
  const double exponent = k[0], max_scoreL = max_score;
  
  splitbit b_compl[MAX_SPLITS][MAX_BINS];
  for (int16 i = 0; i != b.n_splits; i++) {
    for (int16 bin = 0; bin < last_bin; bin++) {
      b_compl[i][bin] = ~b.state[i][bin];
    }
    b_compl[i][last_bin] = b.state[i][last_bin] ^ unset_mask;
  }
  
  int16 a_tips,
  a_and_b, a_and_B, A_and_b, A_and_B,
  a_or_b,  a_or_B,  A_or_b,  A_or_B;
  double ars_ab, ars_aB, ars_Ab, ars_AB,
  min_ars_both, min_ars_either;
  bool enforce_arboreal = arboreal[0];
  
  cost** score = new cost*[most_splits];
  for (int16 i = 0; i < most_splits; i++) score[i] = new cost[most_splits];
  
  for (int16 ai = 0; ai != a.n_splits; ai++) {
    a_tips = 0;
    for (int16 bin = 0; bin != a.n_bins; bin++) {
      a_tips += count_bits(a.state[ai][bin]);
    }
    
    for (int16 bi = 0; bi != b.n_splits; bi++) {
      a_and_b = 0;
      a_and_B = 0;
      a_or_B = 0;
      for (int16 bin = 0; bin != a.n_bins; bin++) {
        a_and_b += count_bits(a.state[ai][bin] & b.state[bi][bin]);
        a_and_B += count_bits(a.state[ai][bin] & b_compl[bi][bin]);
        a_or_B  += count_bits(a.state[ai][bin] | b_compl[bi][bin]);
      }
      A_or_B  = n_tips - a_and_b;
      A_and_b = n_tips - a_or_B;
      A_or_b  = n_tips - a_and_B;
      a_or_b  = a_and_b + a_and_B + A_and_b;
      A_and_B = n_tips - a_or_b;
      
      if (enforce_arboreal && !(
          a_and_b == a_tips ||
            a_and_B == a_tips ||
            A_and_b == n_tips - a_tips ||
            A_and_B == n_tips - a_tips)) {
        
        score[ai][bi] = max_score; /* Prohibit non-arboreal matching */
        
      } else {
        
        ars_ab = (double) a_and_b / (double) a_or_b;
        ars_Ab = (double) A_and_b / (double) A_or_b;
        ars_aB = (double) a_and_B / (double) a_or_B;
        ars_AB = (double) A_and_B / (double) A_or_B;
        
        min_ars_both = (ars_ab < ars_AB) ? ars_ab : ars_AB;
        min_ars_either = (ars_aB < ars_Ab) ? ars_aB : ars_Ab;
        
        /* LAP will look to minimize an integer. max(ars) is between 0 and 1. */
        if (exponent == 1) {
          /* Nye et al. similarity metric */
          score[ai][bi] = cost(max_scoreL - (max_scoreL * 
          ((min_ars_both > min_ars_either) ? 
          min_ars_both : min_ars_either)));
        } else {
          score[ai][bi] = cost(max_scoreL - (max_scoreL * 
            std::pow((min_ars_both > min_ars_either) ? 
            min_ars_both : min_ars_either, exponent)));
        }
      }
    }
    for (int16 bi = b.n_splits; bi < most_splits; bi++) {
      score[ai][bi] = max_score;
    }
  }
  for (int16 ai = a.n_splits; ai < most_splits; ai++) {
    for (int16 bi = 0; bi < most_splits; bi++) {
      score[ai][bi] = max_score;
    }
  }
  
  lap_col *rowsol = new lap_col[most_splits];
  lap_row *colsol = new lap_row[most_splits];
  cost *u = new cost[most_splits], *v = new cost[most_splits];
  
  NumericVector final_score = NumericVector::create(
    (double)((max_score * most_splits) 
               - lap(most_splits, score, rowsol, colsol, u, v))
    / max_score);
  for (int16 i = 0; i < most_splits; i++) delete[] score[i];
  delete[] u; delete[] v; delete[] colsol; delete[] score;
  NumericVector final_matching (most_splits);
  
  for (int16 i = 0; i < most_splits; i++) {
    final_matching[i] = rowsol[i] + 1;
  }
  delete[] rowsol;
  
  List ret = List::create(Named("score") = final_score,
                          _["matching"] = final_matching);
  
  return (ret);
}

// [[Rcpp::export]]
List cpp_mmsi_distance (RawMatrix x, RawMatrix y,
                        IntegerVector nTip) {
  if (x.cols() != y.cols()) {
    throw std::invalid_argument("Input splits must address same number of tips.");
  }
  SplitList a(x), b(y);
  const int16 most_splits = (a.n_splits > b.n_splits) ? a.n_splits : b.n_splits,
    n_tips = nTip[0];
  const cost max_score = BIG;
  const double max_possible = lg2_unrooted[n_tips] - 
    lg2_trees_matching_split((n_tips + 1) / 2, n_tips / 2);
  
  cost** score = new cost*[most_splits];
  for (int16 i = 0; i < most_splits; i++) score[i] = new cost[most_splits];
  
  splitbit different[MAX_BINS];
  int16 n_different, n_same, n_a_only, n_a_and_b;
  double score1, score2;
  
  for (int16 ai = 0; ai != a.n_splits; ai++) {
    for (int16 bi = 0; bi != b.n_splits; bi++) {
      n_different = 0;
      n_a_only = 0;
      n_a_and_b = 0;
      for (int16 bin = 0; bin != a.n_bins; bin++) {
        different[bin] = a.state[ai][bin] ^ b.state[bi][bin];
        n_different += count_bits(different[bin]);
        n_a_only += count_bits(a.state[ai][bin] & different[bin]);
        n_a_and_b += count_bits(a.state[ai][bin] & ~different[bin]);
      }
      n_same = n_tips - n_different;
      
      score1 = lg2_unrooted[n_same] - 
      lg2_trees_matching_split(n_a_and_b, n_same - n_a_and_b);
      
      score2 = lg2_unrooted[n_different] - 
        lg2_trees_matching_split(n_a_only, n_different - n_a_only);
      
      score[ai][bi] = max_score * 
        (1 - ((score1 > score2) ? score1 : score2) / max_possible);
    }
    for (int16 bi = b.n_splits; bi < most_splits; bi++) {
      score[ai][bi] = max_score;
    }
  }
  for (int16 ai = a.n_splits; ai < most_splits; ai++) {
    for (int16 bi = 0; bi < most_splits; bi++) {
      score[ai][bi] = max_score;
    }
  }
  
  lap_col *rowsol = new lap_col[most_splits];
  lap_row *colsol = new lap_row[most_splits];
  cost *u = new cost[most_splits], *v = new cost[most_splits];
  
  NumericVector final_score = NumericVector::create(
    double((max_score * most_splits) - lap(most_splits, score, rowsol, colsol, u, v))
    * max_possible / max_score);
  for (int16 i = 0; i < most_splits; i++) delete[] score[i];
  delete[] u; delete[] v; delete[] colsol; delete[] score;
  NumericVector final_matching (most_splits);
  
  for (int16 i = 0; i < most_splits; i++) {
    final_matching[i] = rowsol[i] + 1;
  }
  delete[] rowsol;
  
  List ret = List::create(Named("score") = final_score,
                          _["matching"] = final_matching);
  
  return (ret);
}

// [[Rcpp::export]]
List cpp_mutual_clustering (RawMatrix x, RawMatrix y,
                            IntegerVector nTip) {
  if (x.cols() != y.cols()) {
    throw std::invalid_argument("Input splits must address same number of tips.");
  }
  SplitList a(x), b(y);
  const int16 most_splits = (a.n_splits > b.n_splits) ? a.n_splits : b.n_splits,
    last_bin = a.n_bins - 1,
    n_tips = nTip[0],
    unset_tips = (n_tips % BIN_SIZE) ? BIN_SIZE - n_tips % BIN_SIZE : 0;
  const cost max_score = BIG;
  const splitbit unset_mask = ALL_ONES >> unset_tips;
  
  splitbit b_compl[MAX_SPLITS][MAX_BINS];
  for (int16 i = 0; i != b.n_splits; i++) {
    for (int16 bin = 0; bin < last_bin; bin++) {
      b_compl[i][bin] = ~b.state[i][bin];
    }
    b_compl[i][last_bin] = b.state[i][last_bin] ^ unset_mask;
  }
  
  cost** score = new cost*[most_splits];
  for (int16 i = 0; i < most_splits; i++) score[i] = new cost[most_splits];
  
  int16 a_and_b, a_and_B, A_and_b, A_and_B,
               na, nA, nb, nB;
  for (int16 ai = 0; ai != a.n_splits; ai++) {
    for (int16 bi = 0; bi != b.n_splits; bi++) {
      a_and_b = 0;
      a_and_B = 0;
      A_and_b = n_tips;
      for (int16 bin = 0; bin != a.n_bins; bin++) {
        a_and_b += count_bits(a.state[ai][bin] & b.state[bi][bin]);
        a_and_B += count_bits(a.state[ai][bin] & b_compl[bi][bin]);
        A_and_b -= count_bits(a.state[ai][bin] | b_compl[bi][bin]);
      }
      A_and_B = n_tips - (a_and_b + a_and_B + A_and_b);
      
      na = a_and_b + a_and_B;
      nA = A_and_b + A_and_B;
      nb = a_and_b + A_and_b;
      nB = a_and_B + A_and_B;
      
      score[ai][bi] = (cost)(max_score * (1L - ((
        ic_element(a_and_b, na, nb, n_tips) +
        ic_element(a_and_B, na, nB, n_tips) +
        ic_element(A_and_b, nA, nb, n_tips) +
        ic_element(A_and_B, nA, nB, n_tips)
      ) / n_tips)));
    }
    for (int16 bi = b.n_splits; bi < most_splits; bi++) {
      score[ai][bi] = max_score;
    }
  }
  for (int16 ai = a.n_splits; ai < most_splits; ai++) {
    for (int16 bi = 0; bi < most_splits; bi++) {
      score[ai][bi] = max_score;
    }
  }
  
  lap_col *rowsol = new lap_col[most_splits];
  lap_row *colsol = new lap_row[most_splits];
  cost *u = new cost[most_splits], *v = new cost[most_splits];
  
  NumericVector final_score = NumericVector::create(
    double((max_score * most_splits) -
      lap(most_splits, score, rowsol, colsol, u, v)) / max_score);
  for (int16 i = 0; i < most_splits; i++) delete[] score[i];
  delete[] colsol; delete[] u; delete[] v; delete[] score;
  
  NumericVector final_matching (most_splits);
  for (int16 i = 0; i != most_splits; i++) {
    final_matching[i] = rowsol[i] + 1;
  }
  delete[] rowsol;
  
  List ret = List::create(Named("score") = final_score,
                          _["matching"] = final_matching);
  
  return (ret);
}

// [[Rcpp::export]]
List cpp_shared_phylo (RawMatrix x, RawMatrix y,
                       IntegerVector nTip) {
  if (x.cols() != y.cols()) {
    throw std::invalid_argument("Input splits must address same number of tips.");
  }
  SplitList a(x), b(y);
  const int16 most_splits = (a.n_splits > b.n_splits) ? a.n_splits : b.n_splits,
              n_tips = nTip[0];
  const cost max_score = BIG;
  const double lg2_unrooted_n = lg2_unrooted[n_tips],
               max_possible = lg2_unrooted_n - 
                 one_overlap((n_tips + 1) / 2, n_tips / 2, n_tips);
  
  int16 in_a[MAX_SPLITS], in_b[MAX_SPLITS];
  for (int16 i = 0; i != a.n_splits; i++) {
    in_a[i] = 0;
    for (int16 bin = 0; bin != a.n_bins; bin++) {
      in_a[i] += count_bits(a.state[i][bin]);
    }
  }
  for (int16 i = 0; i != b.n_splits; i++) {
    in_b[i] = 0;
    for (int16 bin = 0; bin != b.n_bins; bin++) {
      in_b[i] += count_bits(b.state[i][bin]);
    }
  }
  
  cost** score = new cost*[most_splits];
  for (int16 i = 0; i < most_splits; i++) score[i] = new cost[most_splits];
  
  for (int16 ai = 0; ai != a.n_splits; ai++) {
    for (int16 bi = 0; bi != b.n_splits; bi++) {
      score[ai][bi] = max_score * (1 - 
        (spi(a.state[ai], b.state[bi], n_tips, in_a[ai], in_b[bi],
             lg2_unrooted_n, a.n_bins) / max_possible));
    }
    for (int16 bi = b.n_splits; bi < most_splits; bi++) {
      score[ai][bi] = max_score;
    }
  }
  for (int16 ai = a.n_splits; ai < most_splits; ai++) {
    for (int16 bi = 0; bi < most_splits; bi++) {
      score[ai][bi] = max_score;
    }
  }
  
  lap_col *rowsol = new lap_col[most_splits];
  lap_row *colsol = new lap_row[most_splits];
  cost *u = new cost[most_splits], *v = new cost[most_splits];
  
  NumericVector final_score = NumericVector::create(
    (double) ((max_score * most_splits) - lap(most_splits, score, rowsol, colsol, u, v))
    * max_possible / max_score);
  delete[] u; delete[] v; delete[] colsol;
  NumericVector final_matching (most_splits);
  
  
  for (int16 i = 0; i < most_splits; i++) delete[] score[i];
  delete[] score;
  
  for (int16 i = 0; i < most_splits; i++) {
    final_matching[i] = rowsol[i] + 1;
  }
  delete[] rowsol;
  
  List ret = List::create(Named("score") = final_score,
                          _["matching"] = final_matching);
  
  return (ret);
}
