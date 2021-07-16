#include <memory> /* for unique_ptr, make_unique */
#include <cmath>
#include <Rcpp.h>
#include "tree_distances.h"
#include "SplitList.h"

using namespace Rcpp;

// [[Rcpp::export]]
List cpp_robinson_foulds_distance (const RawMatrix x, const RawMatrix y, 
                                   const IntegerVector nTip) {
  if (x.cols() != y.cols()) {
    throw std::invalid_argument("Input splits must address same number of tips.");
  }
  const SplitList a(x), b(y);
  const int16 last_bin = a.n_bins - 1,
              n_tips = nTip[0],
              unset_tips = (n_tips % BIN_SIZE) ? BIN_SIZE - n_tips % BIN_SIZE : 0;
  const splitbit unset_mask = ALL_ONES >> unset_tips;
  cost score = 0;
  
  grf_match matching (a.n_splits);
  for (int16 i = a.n_splits; i--; ) matching[i] = NA_INTEGER;
  
  splitbit b_complement[MAX_SPLITS][MAX_BINS];
  for (int16 i = b.n_splits; i--; ) {
    for (int16 bin = last_bin; bin--; ) {
      b_complement[i][bin] = ~b.state[i][bin];
    }
    b_complement[i][last_bin] = b.state[i][last_bin] ^ unset_mask;
  }
  
  for (int16 ai = a.n_splits; ai--; ) {
    for (int16 bi = b.n_splits; bi--; ) {
    
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
  score = cost(a.n_splits + b.n_splits) - score - score;
  
  return List::create(Named("score") = Rcpp::wrap(score),
                       _["matching"] = Rcpp::wrap(matching));
  
}

// [[Rcpp::export]]
List cpp_robinson_foulds_info (const RawMatrix x, const RawMatrix y, 
                               const IntegerVector nTip) {
  if (x.cols() != y.cols()) {
    throw std::invalid_argument("Input splits must address same number of tips.");
  }
  const SplitList a(x), b(y);
  const int16 last_bin = a.n_bins - 1,
              n_tips = nTip[0],
              unset_tips = (n_tips % BIN_SIZE) ? BIN_SIZE - n_tips % BIN_SIZE : 0;
  const splitbit unset_mask = ALL_ONES >> unset_tips;
  const double lg2_unrooted_n = lg2_unrooted[n_tips];
  double score = 0;
  
  IntegerVector matching (a.n_splits);
  for (int16 i = 0; i != a.n_splits; i++) matching[i] = NA_INTEGER;
  
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
        
        score += lg2_unrooted_n - lg2_rooted[leaves_in_split] -
          lg2_rooted[n_tips - leaves_in_split];

        matching[ai] = bi + 1;
        break; /* Only one match possible per split */
      }
    }
  }
  
  NumericVector final_score = NumericVector::create(score);
  return List::create(Named("score") = final_score,
                      _["matching"] = matching);
  
}

// [[Rcpp::export]]
List cpp_matching_split_distance (const RawMatrix x, const RawMatrix y, 
                                  const IntegerVector nTip) {
  if (x.cols() != y.cols()) {
    throw std::invalid_argument("Input splits must address same number of tips.");
  }
  const SplitList a(x), b(y);
  const int16 most_splits = (a.n_splits > b.n_splits) ? a.n_splits : b.n_splits,
    split_diff = most_splits - 
      ((a.n_splits > b.n_splits) ? b.n_splits : a.n_splits),
    n_tips = nTip[0],
    half_tips = n_tips / 2;
  if (most_splits == 0) {
    return List::create(Named("score") = 0);
  }
  const cost max_score = BIG / most_splits;
  
  cost** score = new cost*[most_splits];
  for (int16 i = most_splits; i--; ) score[i] = new cost[most_splits];
  
  for (int16 ai = 0; ai != a.n_splits; ai++) {
    for (int16 bi = 0; bi != b.n_splits; bi++) {
      score[ai][bi] = 0;
      for (int16 bin = 0; bin != a.n_bins; bin++) {
        score[ai][bi] += count_bits(a.state[ai][bin] ^ b.state[bi][bin]);
      }
      if (score[ai][bi] > half_tips) score[ai][bi] = n_tips - score[ai][bi];
    }
    for (int16 bi = b.n_splits; bi < most_splits; bi++) {
      score[ai][bi] = max_score;
    }
  }
  for (int16 ai = a.n_splits; ai < most_splits; ai++) {
    for (int16 bi = 0; bi != most_splits; bi++) {
      score[ai][bi] = max_score;
    }
  }
  
  lap_col *rowsol = new lap_col[most_splits];
  lap_row *colsol = new lap_row[most_splits];
  cost 
    *u = new cost[most_splits],
    *v = new cost[most_splits]
  ;
  
  NumericVector final_score = NumericVector::create(
    lap(most_splits, score, rowsol, colsol, u, v) - (max_score * split_diff));
  
  for (int16 i = most_splits; i--; ) delete[] score[i];
  delete[] u; delete[] v; delete[] colsol; delete[] score;
  
  IntegerVector final_matching (a.n_splits);
  
  for (int16 i = a.n_splits; i--; ) {
    final_matching[i] = (rowsol[i] < b.n_splits) ? rowsol[i] + 1 : NA_INTEGER;
  }
  
  delete[] rowsol;
  
  return List::create(Named("score") = final_score,
                      _["matching"] = final_matching);
}

// [[Rcpp::export]]
List cpp_jaccard_similarity (const RawMatrix x, const RawMatrix y,
                             const IntegerVector nTip, const NumericVector k,
                             const LogicalVector allowConflict) {
  if (x.cols() != y.cols()) {
    throw std::invalid_argument("Input splits must address same number of tips.");
  }
  const SplitList a(x), b(y);
  const int16
    most_splits = (a.n_splits > b.n_splits) ? a.n_splits : b.n_splits,
    n_tips = nTip[0]
  ;
  const cost max_score = BIG;
  const double exponent = k[0], max_scoreL = max_score;
  
  bool allow_conflict = allowConflict[0];
  
  cost** score = new cost*[most_splits];
  for (int16 i = most_splits; i--; ) score[i] = new cost[most_splits];
  
  for (int16 ai = 0; ai != a.n_splits; ai++) {
    
    const int16 
      na = a.in_split[ai],
      nA = n_tips - na
    ;
    
    for (int16 bi = 0; bi != b.n_splits; bi++) {
      
      // x divides tips into a|A; y divides tips into b|B
      int16 a_and_b = 0;
      for (int16 bin = 0; bin != a.n_bins; bin++) {
        a_and_b += count_bits(a.state[ai][bin] & b.state[bi][bin]);
      }
      
      const int16
        nb = b.in_split[bi],
        nB = n_tips - nb,
        a_and_B = na - a_and_b,
        A_and_b = nb - a_and_b,
        A_and_B = nB - a_and_B
      ;
      
      if (!allow_conflict && !(
            a_and_b == na ||
            a_and_B == na ||
            A_and_b == nA ||
            A_and_B == nA)
          ) {
        
        score[ai][bi] = max_score; /* Prohibited */
        
      } else {
        const int16
          A_or_b  = n_tips - a_and_B,
          a_or_B = n_tips - A_and_b,
          a_or_b = n_tips - A_and_B,
          A_or_B = n_tips - a_and_b
        ;
        const double 
          ars_ab = double(a_and_b) / double(a_or_b),
          ars_Ab = double(A_and_b) / double(A_or_b),
          ars_aB = double(a_and_B) / double(a_or_B),
          ars_AB = double(A_and_B) / double(A_or_B),
          
          min_ars_both = (ars_ab < ars_AB) ? ars_ab : ars_AB,
          min_ars_either = (ars_aB < ars_Ab) ? ars_aB : ars_Ab
        ;
        
        /* LAP will look to minimize an integer. max(ars) is between 0 and 1. */
        if (exponent == 1) {
          /* Nye et al. similarity metric */
          score[ai][bi] = cost(max_scoreL - (max_scoreL * 
            ((min_ars_both > min_ars_either) ? 
            min_ars_both : min_ars_either)));
        } else if (exponent == R_PosInf) {
          score[ai][bi] = cost((min_ars_both == 1 || min_ars_either == 1) ?
                                 0 : max_scoreL);
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
    for (int16 bi = 0; bi != most_splits; bi++) {
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
  for (int16 i = most_splits; i--; ) delete[] score[i];
  delete[] u; delete[] v; delete[] colsol; delete[] score;
  IntegerVector final_matching (a.n_splits);
  
  for (int16 i = a.n_splits; i--; ) {
    final_matching[i] = (rowsol[i] < b.n_splits) ? rowsol[i] + 1 : NA_INTEGER;
  }
  delete[] rowsol;
  
  return List::create(Named("score") = final_score,
                      _["matching"] = final_matching);
  
}

// [[Rcpp::export]]
List cpp_msi_distance (const RawMatrix x, const RawMatrix y,
                        const IntegerVector nTip) {
  if (x.cols() != y.cols()) {
    throw std::invalid_argument("Input splits must address same number of tips.");
  }
  const SplitList a(x), b(y);
  const int16 most_splits = (a.n_splits > b.n_splits) ? a.n_splits : b.n_splits,
              n_tips = nTip[0];
  const cost max_score = BIG;
  const double max_possible = lg2_unrooted[n_tips] - 
    lg2_rooted[int16((n_tips + 1) / 2)] - lg2_rooted[int16(n_tips / 2)];
  
  cost** score = new cost*[most_splits];
  for (int16 i = most_splits; i--; ) score[i] = new cost[most_splits];
  
  splitbit different[MAX_BINS];
  
  for (int16 ai = 0; ai != a.n_splits; ai++) {
    for (int16 bi = 0; bi != b.n_splits; bi++) {
      int16 
        n_different = 0,
        n_a_only = 0,
        n_a_and_b = 0
      ;
      for (int16 bin = 0; bin != a.n_bins; bin++) {
        different[bin] = a.state[ai][bin] ^ b.state[bi][bin];
        n_different += count_bits(different[bin]);
        n_a_only += count_bits(a.state[ai][bin] & different[bin]);
        n_a_and_b += count_bits(a.state[ai][bin] & ~different[bin]);
      }
      const int16 n_same = n_tips - n_different;
      
      score[ai][bi] = max_score - 
        ((max_score / max_possible) *
          mmsi_score(n_same, n_a_and_b, n_different, n_a_only));
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
           lap(most_splits, score, rowsol, colsol, u, v))
    * max_possible / max_score);
  
  for (int16 i = most_splits; i--; ) delete[] score[i];
  delete[] u; delete[] v; delete[] colsol; delete[] score;
  
  IntegerVector final_matching (a.n_splits);
  for (int16 i = a.n_splits; i--; ) {
    final_matching[i] = (rowsol[i] < b.n_splits) ? rowsol[i] + 1 : NA_INTEGER;
  }
  
  delete[] rowsol;
  
  return List::create(Named("score") = final_score,
                      _["matching"] = final_matching);

}

// [[Rcpp::export]]
List cpp_mutual_clustering (const RawMatrix x, const RawMatrix y,
                            const IntegerVector nTip) {
  if (x.cols() != y.cols()) {
    throw std::invalid_argument("Input splits must address same number of tips.");
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
  
  for (int16 ai = 0; ai != a.n_splits; ai++) {
    if (a_match[ai]) continue;
    const int16
      na = a.in_split[ai],
      nA = n_tips - na
    ;
    
    for (int16 bi = 0; bi != b.n_splits; bi++) {
      
      // x divides tips into a|A; y divides tips into b|B
      int16 a_and_b = 0;
      for (int16 bin = 0; bin != a.n_bins; bin++) {
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
    for (int16 bi = b.n_splits; bi < most_splits; bi++) {
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
    for (int16 ai = 0; ai != a.n_splits; ai++) {
      if (a_match[ai]) continue;
      int16 b_pos = 0;
      for (int16 bi = 0; bi != b.n_splits; bi++) {
        if (b_match[bi]) continue;
        score[a_pos][b_pos] = score[ai][bi];
        b_pos++;
      }
      for (int16 bi = lap_dim - a_extra_splits; bi < lap_dim; bi++) {
        score[a_pos][bi] = max_score;
      }
      a_pos++;
    }
    for (int16 ai = lap_dim - b_extra_splits; ai < lap_dim; ai++) {
      for (int16 bi = 0; bi != lap_dim; bi++) {
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
    for (int16 bi = 0; bi != b.n_splits; bi++) {
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
    for (int16 ai = a.n_splits; ai < most_splits; ai++) {
      for (int16 bi = 0; bi != most_splits; bi++) {
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

// [[Rcpp::export]]
List cpp_shared_phylo (const RawMatrix x, const RawMatrix y,
                       const IntegerVector nTip) {
  if (x.cols() != y.cols()) {
    throw std::invalid_argument("Input splits must address same number of tips.");
  }
  const SplitList a(x), b(y);
  const int16
    most_splits = (a.n_splits > b.n_splits) ? a.n_splits : b.n_splits,
    n_tips = nTip[0]
  ;
  const cost max_score = BIG;
  const double
    lg2_unrooted_n = lg2_unrooted[n_tips],
    best_overlap = one_overlap((n_tips + 1) / 2, n_tips / 2, n_tips),
    max_possible = lg2_unrooted_n - best_overlap
  ;
  
  // a and b are "clades" separating an "ingroup" [1] from an "outgroup" [0].
  // In/out direction [i.e. 1/0 bit] is arbitrary.
  cost** score = new cost*[most_splits];
  for (int16 i = most_splits; i--; ) score[i] = new cost[most_splits];
  
  for (int16 ai = a.n_splits; ai--; ) {
    for (int16 bi = b.n_splits; bi--; ) {
      const double spi_over = spi_overlap(a.state[ai], b.state[bi], n_tips,
                                          a.in_split[ai], b.in_split[bi],
                                          a.n_bins);
      
      score[ai][bi] = spi_over ?
        (spi_over - best_overlap) * (max_score / max_possible) :
        max_score;
        
    }
    for (int16 bi = b.n_splits; bi < most_splits; ++bi) {
      score[ai][bi] = max_score;
    }
  }
  for (int16 ai = a.n_splits; ai < most_splits; ++ai) {
    for (int16 bi = 0; bi != most_splits; ++bi) {
      score[ai][bi] = max_score;
    }
  }
  
  lap_col *rowsol = new lap_col[most_splits];
  lap_row *colsol = new lap_row[most_splits];
  cost *u = new cost[most_splits], *v = new cost[most_splits];
  
  NumericVector final_score = NumericVector::create(
    double((max_score * most_splits) - 
      lap(most_splits, score, rowsol, colsol, u, v))
    * (max_possible / max_score));
  
  delete[] u; delete[] v; delete[] colsol;
  
  IntegerVector final_matching (a.n_splits);
  
  for (int16 i = most_splits; i--; ) delete[] score[i];
  delete[] score;
  
  for (int16 i = a.n_splits; i--; ) {
    final_matching[i] = (rowsol[i] < b.n_splits) ? rowsol[i] + 1 : NA_INTEGER;
  }
  
  delete[] rowsol;
  
  return List::create(Named("score") = final_score,
                      _["matching"] = final_matching);

}
 