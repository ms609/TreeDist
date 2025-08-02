#include <TreeTools/SplitList.h>
#include <cmath>
#include <memory> /* for unique_ptr, make_unique */
#include <Rcpp/Lightest>
#include "tree_distances.h"

using namespace Rcpp;
using TreeTools::SplitList;
using TreeTools::count_bits;

namespace TreeDist {

  void check_ntip(const double n) {
    if (n > static_cast<double>(std::numeric_limits<int16>::max())) {
      Rcpp::stop("This many tips are not (yet) supported.");
    }
  }

  inline void add_ic_element(double& ic_sum, int16 nkK, int16 nk, int16 nK,
                             int16 n_tips) noexcept {
    /* 
     * See equation 16 in Meila 2007. I denote k' as K.
     * nkK is converted to pkK in the calling function, when the sum of all
     * elements is divided by n.
     */
    if (nkK && nk && nK) {
      if (nkK == nk && nkK == nK && nkK << 1 == n_tips) {
        ic_sum += nkK;
      } else {
        const int32 numerator = nkK * n_tips;
        const int32 denominator = nk * nK;
        if (numerator != denominator) {
          ic_sum += nkK * (lg2[numerator] - lg2[denominator]);
        }
      }
    }
  }

}



// [[Rcpp::export]]
List cpp_robinson_foulds_distance(const RawMatrix x, const RawMatrix y, 
                                  const IntegerVector nTip) {
  if (x.cols() != y.cols()) {
    Rcpp::stop("Input splits must address same number of tips.");
  }
  TreeDist::check_ntip(nTip[0]);
  
  const SplitList a(x), b(y);
  const int16 last_bin = a.n_bins - 1;
  const int16 n_tips = int16(nTip[0]);
  const int16 unset_tips = (n_tips % SL_BIN_SIZE) ?
    SL_BIN_SIZE - n_tips % SL_BIN_SIZE : 0;
  const splitbit unset_mask = ALL_ONES >> unset_tips;
  cost score = 0;
  
  grf_match matching(a.n_splits, NA_INTEGER);
  
  splitbit b_complement[SL_MAX_SPLITS][SL_MAX_BINS];
  for (int16 i = b.n_splits; i--; ) {
    for (int16 bin = last_bin; bin--; ) {
      b_complement[i][bin] = ~b.state[i][bin];
    }
    b_complement[i][last_bin] = b.state[i][last_bin] ^ unset_mask;
  }
  
  for (int16 ai = a.n_splits; ai--; ) {
    for (int16 bi = b.n_splits; bi--; ) {
      
      bool all_match = true;
      bool all_complement = true;
      
      for (int16 bin = 0; bin < a.n_bins; ++bin) {
        if ((a.state[ai][bin] != b.state[bi][bin])) {
          all_match = false;
          break;
        }
      }
      if (!all_match) {
        for (int16 bin = 0; bin < a.n_bins; ++bin) {
          if (a.state[ai][bin] != b_complement[bi][bin]) {
            all_complement = false;
            break;
          }
        }
      }
      if (all_match || all_complement) {
        ++score;
        matching[ai] = bi + 1;
        break; // Only one match possible per split
      }
    }
  }
  score = cost(a.n_splits + b.n_splits) - score - score;
  
  return List::create(Named("score") = Rcpp::wrap(score),
                      _["matching"] = Rcpp::wrap(matching));
  
}

// [[Rcpp::export]]
List cpp_robinson_foulds_info(const RawMatrix x, const RawMatrix y, 
                              const IntegerVector nTip) {
  if (x.cols() != y.cols()) {
    Rcpp::stop("Input splits must address same number of tips.");
  }
  TreeDist::check_ntip(nTip[0]);
  
  const SplitList a(x), b(y);
  
  const int16 last_bin = a.n_bins - 1;
  const int16 n_tips = int16(nTip[0]);
  const int16 unset_tips = (n_tips % SL_BIN_SIZE) ?
    SL_BIN_SIZE - n_tips % SL_BIN_SIZE : 0;
  
  const splitbit unset_mask = ALL_ONES >> unset_tips;
  const double lg2_unrooted_n = lg2_unrooted[n_tips];
  double score = 0;
  
  grf_match matching(a.n_splits, NA_INTEGER);
  
  /* Dynamic allocation 20% faster for 105 tips, but VLA not permitted in C11 */
  splitbit b_complement[SL_MAX_SPLITS][SL_MAX_BINS]; 
  for (int16 i = 0; i < b.n_splits; i++) {
    for (int16 bin = 0; bin < last_bin; ++bin) {
      b_complement[i][bin] = ~b.state[i][bin];
    }
    b_complement[i][last_bin] = b.state[i][last_bin] ^ unset_mask;
  }
  
  for (int16 ai = 0; ai < a.n_splits; ++ai) {
    for (int16 bi = 0; bi < b.n_splits; ++bi) {
      
      bool all_match = true, all_complement = true;
      
      for (int16 bin = 0; bin < a.n_bins; ++bin) {
        if ((a.state[ai][bin] != b.state[bi][bin])) {
          all_match = false;
          break;
        }
      }
      if (!all_match) {
        for (int16 bin = 0; bin < a.n_bins; ++bin) {
          if ((a.state[ai][bin] != b_complement[bi][bin])) {
            all_complement = false;
            break;
          }
        }
      }
      if (all_match || all_complement) {
        int16 leaves_in_split = 0;
        for (int16 bin = 0; bin < a.n_bins; ++bin) {
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
List cpp_matching_split_distance(const RawMatrix x, const RawMatrix y, 
                                 const IntegerVector nTip) {
  if (x.cols() != y.cols()) {
    Rcpp::stop("Input splits must address same number of tips.");
  }
  TreeDist::check_ntip(nTip[0]);
  
  const SplitList a(x), b(y);
  const int16 most_splits = (a.n_splits > b.n_splits) ? a.n_splits : b.n_splits;
  const int16 split_diff = most_splits - 
    ((a.n_splits > b.n_splits) ? b.n_splits : a.n_splits);
  const int16 n_tips = int16(nTip[0]);
  const int16 half_tips = n_tips / 2;
  if (most_splits == 0) {
    return List::create(Named("score") = 0);
  }
  const cost max_score = BIG / most_splits;
  cost_matrix score(most_splits);
  
  for (int16 ai = 0; ai < a.n_splits; ++ai) {
    for (int16 bi = 0; bi < b.n_splits; ++bi) {
      splitbit total = 0;
      for (int16 bin = 0; bin < a.n_bins; ++bin) { // hot hot hot 11.3% -> 20.4%
         total += count_bits(a.state[ai][bin] ^ b.state[bi][bin]); // hot hot 8.0% → 15%
      }
      score(ai, bi) = total;
    }
  }
  
  // More performat to cleave this out, to keep the previous loop simpler
  for (int16 ai = 0; ai < a.n_splits; ++ai) {
    for (int16 bi = 0; bi < b.n_splits; ++bi) {
      if (score(ai, bi) > half_tips) {
        score(ai, bi) = n_tips - score(ai, bi);
      }
    }
    score.padRowAfterCol(ai, b.n_splits, max_score);
  }
  
  score.padAfterRow(a.n_splits, max_score);
  
  std::vector<lap_col> rowsol(most_splits);
  std::vector<lap_row> colsol(most_splits);
  
  NumericVector final_score = NumericVector::create(
    lap(most_splits, score, rowsol, colsol) - (max_score * split_diff));
  
  IntegerVector final_matching(a.n_splits);
  
  for (int16 i = 0; i < a.n_splits; ++i) {
    final_matching[i] = (rowsol[i] < b.n_splits) ? rowsol[i] + 1 : NA_INTEGER;
  }
  
  return List::create(Named("score") = final_score,
                      _["matching"] = final_matching);
}

// [[Rcpp::export]]
List cpp_jaccard_similarity(const RawMatrix x, const RawMatrix y,
                            const IntegerVector nTip, const NumericVector k,
                            const LogicalVector allowConflict) {
  if (x.cols() != y.cols()) {
    Rcpp::stop("Input splits must address same number of tips.");
  }
  TreeDist::check_ntip(nTip[0]);
  
  const SplitList a(x), b(y);
  const int16 most_splits = (a.n_splits > b.n_splits) ? a.n_splits : b.n_splits;
  const int16 n_tips = int16(nTip[0]);
  
  constexpr cost max_score = BIG;
  constexpr double max_scoreL = max_score;
  const double exponent = k[0];
  
  bool allow_conflict = allowConflict[0];
  
  cost_matrix score(most_splits);
  
  for (int16 ai = 0; ai < a.n_splits; ++ai) {
    
    const int16 na = a.in_split[ai];
    const int16 nA = n_tips - na;
    
    for (int16 bi = 0; bi < b.n_splits; ++bi) {
      
      // x divides tips into a|A; y divides tips into b|B
      int16 a_and_b = 0;
      for (int16 bin = 0; bin < a.n_bins; ++bin) {
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
        
        score(ai, bi) = max_score; /* Prohibited */
        
      } else {
        const int16
          A_or_b  = n_tips - a_and_B,
          a_or_B = n_tips - A_and_b,
          a_or_b = n_tips - A_and_B,
          A_or_B = n_tips - a_and_b
        ;
        const double 
          ars_ab = static_cast<double>(a_and_b) / static_cast<double>(a_or_b),
          ars_Ab = static_cast<double>(A_and_b) / static_cast<double>(A_or_b),
          ars_aB = static_cast<double>(a_and_B) / static_cast<double>(a_or_B),
          ars_AB = static_cast<double>(A_and_B) / static_cast<double>(A_or_B),
          
          min_ars_both = (ars_ab < ars_AB) ? ars_ab : ars_AB,
          min_ars_either = (ars_aB < ars_Ab) ? ars_aB : ars_Ab
        ;
        
        /* LAP will look to minimize an integer. max(ars) is between 0 and 1. */
        if (exponent == 1) {
          /* Nye et al. similarity metric */
          score(ai, bi) = cost(max_scoreL - (max_scoreL * 
            ((min_ars_both > min_ars_either) ? 
            min_ars_both : min_ars_either)));
        } else if (exponent == R_PosInf) {
          score(ai, bi) = cost((min_ars_both == 1 || min_ars_either == 1) ?
                                 0 : max_scoreL);
        } else {
          score(ai, bi) = cost(max_scoreL - (max_scoreL * 
            std::pow((min_ars_both > min_ars_either) ? 
            min_ars_both : min_ars_either, exponent)));
        }
      }
    }
    score.padRowAfterCol(ai, b.n_splits, max_score);
  }
  score.padAfterRow(a.n_splits, max_score);
  
  std::vector<lap_col> rowsol(most_splits);
  std::vector<lap_row> colsol(most_splits);
  
  NumericVector final_score = NumericVector::create(
    (double)((max_score * most_splits) - 
      lap(most_splits, score, rowsol, colsol))
    / max_score);
  IntegerVector final_matching(a.n_splits);
  
  for (int16 i = 0; i < a.n_splits; ++i) {
    final_matching[i] = (rowsol[i] < b.n_splits) ? rowsol[i] + 1 : NA_INTEGER;
  }
  
  return List::create(Named("score") = final_score,
                      _["matching"] = final_matching);
  
}

// [[Rcpp::export]]
List cpp_msi_distance(const RawMatrix x, const RawMatrix y,
                      const IntegerVector nTip) {
  if (x.cols() != y.cols()) {
    Rcpp::stop("Input splits must address same number of tips.");
  }
  TreeDist::check_ntip(nTip[0]);
  
  const SplitList a(x), b(y);
  const int16 most_splits = (a.n_splits > b.n_splits) ? a.n_splits : b.n_splits;
  const int16 n_tips = int16(nTip[0]);
  constexpr cost max_score = BIG;
  const double max_possible = lg2_unrooted[n_tips] - 
    lg2_rooted[int16((n_tips + 1) / 2)] - lg2_rooted[int16(n_tips / 2)];
  const double score_over_possible = static_cast<double>(max_score) / max_possible;
  const double possible_over_score = max_possible / max_score;
  
  cost_matrix score(most_splits);
  
  splitbit different[SL_MAX_BINS];
  
  for (int16 ai = 0; ai < a.n_splits; ++ai) {
    for (int16 bi = 0; bi < b.n_splits; ++bi) {
      int16 
        n_different = 0,
        n_a_only = 0,
        n_a_and_b = 0
      ;
      for (int16 bin = 0; bin < a.n_bins; ++bin) {
        different[bin] = a.state[ai][bin] ^ b.state[bi][bin];
        n_different += count_bits(different[bin]);
        n_a_only += count_bits(a.state[ai][bin] & different[bin]);
        n_a_and_b += count_bits(a.state[ai][bin] & ~different[bin]);
      }
      const int16 n_same = n_tips - n_different;
      
      score(ai, bi) = cost(max_score - 
        (score_over_possible *
          TreeDist::mmsi_score(n_same, n_a_and_b, n_different, n_a_only)));
    }
    score.padRowAfterCol(ai, b.n_splits, max_score);
  }
  score.padAfterRow(a.n_splits, max_score);
  
  std::vector<lap_col> rowsol(most_splits);
  std::vector<lap_row> colsol(most_splits);
  
  NumericVector final_score = NumericVector::create(
    static_cast<double>((max_score * most_splits) - 
                        lap(most_splits, score, rowsol, colsol))
    * possible_over_score
  );
  
  IntegerVector final_matching(a.n_splits);
  for (int16 i = 0; i < a.n_splits; ++i) {
    final_matching[i] = (rowsol[i] < b.n_splits) ? rowsol[i] + 1 : NA_INTEGER;
  }
  
  return List::create(Named("score") = final_score,
                      _["matching"] = final_matching);

}

// [[Rcpp::export]]
List cpp_mutual_clustering(const RawMatrix x, const RawMatrix y,
                           const IntegerVector nTip) {
  if (x.cols() != y.cols()) {
    Rcpp::stop("Input splits must address same number of tips.");
  }
  TreeDist::check_ntip(nTip[0]);
  
  const SplitList a(x);
  const SplitList b(y);
  const bool a_has_more_splits = (a.n_splits > b.n_splits);
  const int16 most_splits = a_has_more_splits ? a.n_splits : b.n_splits;
  const int16 a_extra_splits = a_has_more_splits ? most_splits - b.n_splits : 0;
  const int16 b_extra_splits = a_has_more_splits ? 0 : most_splits - a.n_splits;
  const int16 n_tips = int16(nTip[0]);
  const double n_tips_reciprocal = 1.0 / n_tips;
  
  if (most_splits == 0 || n_tips == 0) {
    return List::create(Named("score") = 0,
                        _["matching"] = IntegerVector(0));
  }
  constexpr cost max_score = BIG;
  const double max_over_tips = static_cast<double>(max_score) * n_tips_reciprocal;
  cost_matrix score(most_splits);
  
  double exact_match_score = 0;
  int16 exact_matches = 0;
  // NumericVector zero-initializes [so does make_unique]
  // match will have one added to it so numbering follows R; hence 0 = UNMATCHED
  IntegerVector a_match(a.n_splits);
  std::unique_ptr<int16[]> b_match = std::make_unique<int16[]>(b.n_splits);
  
  for (int16 ai = 0; ai < a.n_splits; ++ai) {
    if (a_match[ai]) continue;
    const int16 na = a.in_split[ai];
    const int16 nA = n_tips - na;
    
    for (int16 bi = 0; bi < b.n_splits; ++bi) {
      
      // x divides tips into a|A; y divides tips into b|B
      int16 a_and_b = 0;
      for (int16 bin = 0; bin < a.n_bins; ++bin) {
        a_and_b += count_bits(a.state[ai][bin] & b.state[bi][bin]);
      }
      
      const int16 nb = b.in_split[bi];
      const int16 nB = n_tips - nb;
      const int16 a_and_B = na - a_and_b;
      const int16 A_and_b = nb - a_and_b;
      const int16 A_and_B = nA - A_and_b;
      
      if ((!a_and_B && !A_and_b) ||
          (!a_and_b && !A_and_B)) {
        exact_match_score += TreeDist::ic_matching(na, nA, n_tips);
        ++exact_matches;
        a_match[ai] = bi + 1;
        b_match[bi] = ai + 1;
        break;
      } else if (a_and_b == A_and_b &&
          a_and_b == a_and_B &&
          a_and_b == A_and_B) {
        score(ai, bi) = max_score; // Don't risk rounding error
      } else {
        double ic_sum = 0;
        TreeDist::add_ic_element(ic_sum, a_and_b, na, nb, n_tips);
        TreeDist::add_ic_element(ic_sum, a_and_B, na, nB, n_tips);
        TreeDist::add_ic_element(ic_sum, A_and_b, nA, nb, n_tips);
        TreeDist::add_ic_element(ic_sum, A_and_B, nA, nB, n_tips);
        
        // Division by n_tips converts n(A&B) to P(A&B) for each ic_element
        score(ai, bi) = max_score - static_cast<cost>(ic_sum * max_over_tips);
      }
    }
    score.padRowAfterCol(ai, b.n_splits, max_score);
  }
  if (exact_matches == b.n_splits || exact_matches == a.n_splits) {
    return List::create(
      Named("score") = NumericVector::create(exact_match_score * n_tips_reciprocal),
      _["matching"] = a_match);
  }
  
  
  const int16 lap_dim = most_splits - exact_matches;
  std::vector<lap_col> rowsol(lap_dim);
  std::vector<lap_row> colsol(lap_dim);
  cost_matrix small_score(lap_dim);
  
  if (exact_matches) {
    int16 a_pos = 0;
    for (int16 ai = 0; ai < a.n_splits; ++ai) {
      if (a_match[ai]) continue;
      int16 b_pos = 0;
      for (int16 bi = 0; bi < b.n_splits; ++bi) {
        if (b_match[bi]) continue;
        small_score(a_pos, b_pos) = score(ai, bi);
        b_pos++;
      }
      small_score.padRowAfterCol(a_pos, lap_dim - a_extra_splits, max_score);
      a_pos++;
    }
    small_score.padAfterRow(lap_dim - b_extra_splits, max_score);
    
    const double lap_score = static_cast<double>(
      (max_score * lap_dim) - lap(lap_dim, small_score, rowsol, colsol)
      ) / max_score;
    NumericVector final_score = 
      NumericVector::create(lap_score + (exact_match_score * n_tips_reciprocal));
    
    std::unique_ptr<int16[]> lap_decode = std::make_unique<int16[]>(lap_dim);
    int16 fuzzy_match = 0;
    for (int16 bi = 0; bi < b.n_splits; ++bi) {
      if (!b_match[bi]) {
        assert(fuzzy_match < lap_dim);
        lap_decode[fuzzy_match++] = bi + 1;
      }
    }
    
    fuzzy_match = 0;
    IntegerVector final_matching(a.n_splits);
    for (int16 i = 0; i < a.n_splits; ++i) {
      if (a_match[i]) {
        // Rcout << "a" << (1+i) << " exactly matches b" << a_match[i]<< "\n";
        final_matching[i] = a_match[i];
      } else {
        assert(fuzzy_match < lap_dim);
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
    
    return List::create(Named("score") = final_score,
                        _["matching"] = final_matching);
  } else {
    for (int16 ai = a.n_splits; ai < most_splits; ++ai) {
      for (int16 bi = 0; bi < most_splits; ++bi) {
        score(ai, bi) = max_score;
      }
    }
    
    const double lap_score = static_cast<double>(
        (max_score * lap_dim) - lap(lap_dim, score, rowsol, colsol)
      ) / max_score;
    NumericVector final_score = NumericVector::create(lap_score);
    
    IntegerVector final_matching(a.n_splits);
    for (int16 i = 0; i < a.n_splits; ++i) {
      final_matching[i] = (rowsol[i] < b.n_splits) ? rowsol[i] + 1 : NA_INTEGER;
    }
    
    return List::create(Named("score") = final_score,
                        _["matching"] = final_matching);
  }
}

// [[Rcpp::export]]
List cpp_shared_phylo (const RawMatrix x, const RawMatrix y,
                       const IntegerVector nTip) {
  if (x.cols() != y.cols()) {
    Rcpp::stop("Input splits must address same number of tips.");
  }
  TreeDist::check_ntip(nTip[0]);
  
  const SplitList a(x), b(y);
  const int16 most_splits = (a.n_splits > b.n_splits) ? a.n_splits : b.n_splits;
  const int16 n_tips = int16(nTip[0]);
  const int16 overlap_a = int16(n_tips + 1) / 2; // avoids promotion to int
  
  constexpr cost max_score = BIG;
  const double lg2_unrooted_n = lg2_unrooted[n_tips];
  const double best_overlap = TreeDist::one_overlap(overlap_a, n_tips / 2, n_tips);
  const double max_possible = lg2_unrooted_n - best_overlap;
  const double score_over_possible = max_score / max_possible;
  const double possible_over_score = max_possible / max_score;
  
  // a and b are "clades" separating an "ingroup" [1] from an "outgroup" [0].
  // In/out direction [i.e. 1/0 bit] is arbitrary.
  cost_matrix score(most_splits);
  
  for (int16 ai = a.n_splits; ai--; ) {
    for (int16 bi = b.n_splits; bi--; ) {
      const double spi_over = TreeDist::spi_overlap(
        a.state[ai], b.state[bi], n_tips, a.in_split[ai], b.in_split[bi],
                                                                    a.n_bins);
      
      score(ai, bi) = spi_over == 0 ? max_score :
        cost((spi_over - best_overlap) * score_over_possible);
    }
    score.padRowAfterCol(ai, b.n_splits, max_score);
  }
  score.padAfterRow(a.n_splits, max_score);
  
  std::vector<lap_col> rowsol(most_splits);
  std::vector<lap_row> colsol(most_splits);
  
  NumericVector final_score = NumericVector::create(
    static_cast<double>(
      (max_score * most_splits) - lap(most_splits, score, rowsol, colsol)) *
        possible_over_score
    );
  
  IntegerVector final_matching(a.n_splits);
  
  for (int16 i = 0; i < a.n_splits; ++i) {
    final_matching[i] = (rowsol[i] < b.n_splits) ? rowsol[i] + 1 : NA_INTEGER;
  }
  
  return List::create(Named("score") = final_score,
                      _["matching"] = final_matching);

}
