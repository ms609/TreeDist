#include <Rcpp.h>
using namespace Rcpp;
#include <stdint.h>
#include <math.h> /* for pow(), log() */
#include "SplitList.h"
#include "lap.h"

uint32_t bitcounts[65536]; // the bytes representing bit count of each number 0-65535
__attribute__((constructor))
  void initialize_bitcounts() {
    for (int32_t i = 0; i < 65536; i++) {
      int32_t n_bits = 0;
      for (int j = 0; j < 16; j++) {
        if ((i & powers_of_two[j])) ++n_bits;
      }
      bitcounts[i] = n_bits;
    }
  }

int count_bits_32 (uint32_t x) {
  return bitcounts[x & right16bits] + bitcounts[x >> 16];
}

double lg2_double_factorial[6398]; /* Only offer support for up to 3200 tips */
double lg2_rooted[3201];
double lg2_unrooted[3201];
__attribute__((constructor))
  void initialize_ldf() {
    for (int i = 0; i < 3; i++) {
      lg2_double_factorial[i] = 0;
      lg2_rooted[i] = 0;
      lg2_unrooted[i] = 0;
    }
    for (int i = 2; i < 6398; i++) {
      lg2_double_factorial[i] = lg2_double_factorial[i - 2] + log2(i);
    }
    for (int i = 3; i < 3201; i++) {
      lg2_unrooted[i] = lg2_double_factorial[i + i - 5];
      lg2_rooted[i] = lg2_double_factorial[i + i - 3];
    }
  }

// [[Rcpp::export]]
List cpp_robinson_foulds_distance (NumericMatrix x, NumericMatrix y, 
                                  NumericVector nTip) {
  if (x.cols() != y.cols()) {
    throw std::invalid_argument("Input splits must address same number of tips.");
  }
  SplitList a(x), b(y);
  const int max_splits = (a.n_splits > b.n_splits) ? a.n_splits : b.n_splits,
            last_bin = a.n_bins - 1,
            n_tips = nTip[0],
            unset_tips = (n_tips % 32) ? 32 - n_tips % 32 : 32;
  const uint32_t unset_mask = ~0U >> unset_tips;
  
  int score = 0;
  NumericVector matching (max_splits);
  for (int i = 0; i < max_splits; i++) matching[i] = NA_REAL;
  
  uint32_t b_complement[b.n_splits][b.n_bins];
  for (int i = 0; i < b.n_splits; i++) {
    for (int bin = 0; bin < last_bin; bin++) {
        b_complement[i][bin] = ~b.state[i][bin];
    }
    b_complement[i][last_bin] = b.state[i][last_bin] ^ unset_mask;
  }
  
  for (int ai = 0; ai < a.n_splits; ai++) {
    for (int bi = 0; bi < b.n_splits; bi++) {
      bool all_match = true, all_mismatch = true;
      for (int bin = 0; bin < a.n_bins; bin++) {
        if ((a.state[ai][bin] != b.state[bi][bin])) {
          all_match = false;
          break;
        }
      }
      if (!all_match) {
        for (int bin = 0; bin < a.n_bins; bin++) {
          if ((a.state[ai][bin] != b_complement[bi][bin])) {
            all_mismatch = false;
            break;
          }
        }
      }
      if (all_match || all_mismatch) {
        ++score;
        matching[ai] = bi;
        break; /* Only one match possible per split */
      }
    }
  }
  score = a.n_splits + b.n_splits - score - score;
  
  NumericVector final_score = NumericVector::create(score);
  
  List ret = List::create(Named("score") = final_score,
                          _["matching"] = matching);
  
  return (ret);
}

// [[Rcpp::export]]
List cpp_matching_split_distance (NumericMatrix x, NumericMatrix y, 
                                  NumericVector nTip) {
  if (x.cols() != y.cols()) {
    throw std::invalid_argument("Input splits must address same number of tips.");
  }
  SplitList a(x), b(y);
  const int max_splits = (a.n_splits > b.n_splits) ? a.n_splits : b.n_splits,
            split_diff = max_splits - 
                         ((a.n_splits > b.n_splits) ? b.n_splits : a.n_splits),
            n_tips = nTip[0],
            half_tips = n_tips / 2;
  
  int** score = new int*[max_splits];
  for (int i = 0; i < max_splits; i++) score[i] = new int[max_splits];
  
  /*Rcout << "Working over " << a.n_splits << " (" << a.n_splits << ", " << x.rows() 
        << ") and " << b.n_splits << " (" << b.n_splits << ", " << y.rows() 
        << ") splits.\n\n";*/
  
  for (int ai = 0; ai < a.n_splits; ai++) {
    for (int bi = 0; bi < b.n_splits; bi++) {
      score[ai][bi] = 0;
      for (int bin = 0; bin < a.n_bins; bin++) {
        score[ai][bi] += count_bits_32(a.state[ai][bin] ^ 
                                       b.state[bi][bin]);
        /*Rcout << "- x = " << ai << ", y = " << bi << ", bin " << bin << ": "
              << a.state[ai][bin] << " ^ " << b.state[bi][bin] << " = " 
              << score[ai][bi] << " (" << n_tips << " tips).\n";*/
      }
      if (score[ai][bi] > half_tips) score[ai][bi] = n_tips - score[ai][bi];
    }
    for (int bi = b.n_splits; bi < max_splits; bi++) {
      score[ai][bi] = BIG;
    }
  }
  for (int ai = a.n_splits; ai < max_splits; ai++) {
    for (int bi = 0; bi < max_splits; bi++) {
      score[ai][bi] = BIG;
    }
  }
  
  lap_col *rowsol = new lap_col[max_splits];
  lap_row *colsol = new lap_row[max_splits];
  cost *u = new cost[max_splits], 
       *v = new cost[max_splits];
  
  NumericVector final_score = NumericVector::create(
    lap(max_splits, score, rowsol, colsol, u, v) - (BIG * split_diff)),
    final_matching (max_splits);
  
  for (int i = 0; i < max_splits; i++) {
    final_matching[i] = rowsol[i] + 1;
  }
  
  List ret = List::create(Named("score") = final_score,
                          _["matching"] = final_matching);
  
  return (ret);
}

// [[Rcpp::export]]
List cpp_jaccard_distance (NumericMatrix x, NumericMatrix y,
                           NumericVector nTip, NumericVector k,
                           LogicalVector arboreal) {
  if (x.cols() != y.cols()) {
    throw std::invalid_argument("Input splits must address same number of tips.");
  }
  SplitList a(x), b(y);
  const int max_splits = (a.n_splits > b.n_splits) ? a.n_splits : b.n_splits,
    last_bin = a.n_bins - 1,
    n_tips = nTip[0],
    unset_tips = (n_tips % 32) ? 32 - n_tips % 32 : 32;
  const uint32_t unset_mask = ~0U >> unset_tips;
  const double exponent = k[0];
  
  uint32_t b_compl[b.n_splits][b.n_bins];
  for (int i = 0; i < b.n_splits; i++) {
    for (int bin = 0; bin < last_bin; bin++) {
      b_compl[i][bin] = ~b.state[i][bin];
    }
    b_compl[i][last_bin] = b.state[i][last_bin] ^ unset_mask;
  }
  
  int a_and_b, a_and_B, A_and_b, A_and_B,
    a_or_b,  a_or_B,  A_or_b,  A_or_B;
  double ars_ab, ars_aB, ars_Ab, ars_AB,
    min_ars_both, min_ars_either;
  bool enforce_arboreal = arboreal[0];
  
  int** score = new int*[max_splits];
  for (int i = 0; i < max_splits; i++) score[i] = new int[max_splits];
  
  for (int ai = 0; ai < a.n_splits; ai++) {
    for (int bi = 0; bi < b.n_splits; bi++) {
      a_and_b = 0;
      a_and_B = 0;
      a_or_B = 0;
      for (int bin = 0; bin < a.n_bins; bin++) {
        a_and_b += count_bits_32(a.state[ai][bin] & b.state[bi][bin]);
        a_and_B += count_bits_32(a.state[ai][bin] & b_compl[bi][bin]);
        a_or_B  += count_bits_32(a.state[ai][bin] | b_compl[bi][bin]);
      }
      A_or_B  = n_tips - a_and_b;
      A_and_b = n_tips - a_or_B;
      A_or_b  = n_tips - a_and_B;
      a_or_b  = a_and_b + a_and_B + A_and_b;
      A_and_B = n_tips - a_or_b;
      
      if (enforce_arboreal && !(
          a_and_b == n_tips ||
            a_and_B == n_tips ||
            A_and_b == n_tips ||
            A_and_B == n_tips)) {
        
        score[ai][bi] = BIG; /* Prohibit non-arboreal matching */
        
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
          score[ai][bi] = (int) BIGL - (BIGL * 
          ((min_ars_both > min_ars_either) ? 
          min_ars_both : min_ars_either));
        } else {
          /*Rcout << "Score: " << ((min_ars_both > min_ars_either) ? 
           min_ars_both : min_ars_either)
           << " ^ " << exponent << " = " << pow((min_ars_both > min_ars_either) ? 
           min_ars_both : min_ars_either, exponent) 
           << ", BIG - BIG*score = " <<( (int) BIGL - (BIGL * 
           pow((min_ars_both > min_ars_either) ? 
           min_ars_both : min_ars_either, exponent))) << ".\n";*/
          score[ai][bi] = (int) BIGL - (BIGL * 
            pow((min_ars_both > min_ars_either) ? 
            min_ars_both : min_ars_either, exponent));
        }
      }
    }
    for (int bi = b.n_splits; bi < max_splits; bi++) {
      score[ai][bi] = BIG;
    }
  }
  for (int ai = a.n_splits; ai < max_splits; ai++) {
    for (int bi = 0; bi < max_splits; bi++) {
      score[ai][bi] = BIG;
    }
  }
  
  lap_col *rowsol = new lap_col[max_splits];
  lap_row *colsol = new lap_row[max_splits];
  cost *u = new cost[max_splits], *v = new cost[max_splits];
  
  NumericVector final_score = NumericVector::create(
    (double)((BIG * max_splits) - lap(max_splits, score, rowsol, colsol, u, v))
    / BIGL),
    final_matching (max_splits);
  
  for (int i = 0; i < max_splits; i++) {
    final_matching[i] = rowsol[i] + 1;
  }
  
  List ret = List::create(Named("score") = final_score,
                          _["matching"] = final_matching);
  
  return (ret);
}

double lg2_trees_matching_split (int a, int b) {
  if (a == 0) return (lg2_unrooted[b]);
  if (b == 0) return (lg2_unrooted[a]);
  return(lg2_rooted[a] + lg2_rooted[b]);
}

// [[Rcpp::export]]
List cpp_mmsi_distance (NumericMatrix x, NumericMatrix y,
                           NumericVector nTip) {
  if (x.cols() != y.cols()) {
    throw std::invalid_argument("Input splits must address same number of tips.");
  }
  SplitList a(x), b(y);
  const int max_splits = (a.n_splits > b.n_splits) ? a.n_splits : b.n_splits,
    n_tips = nTip[0];
  const double max_score = lg2_unrooted[n_tips] - 
    lg2_trees_matching_split((n_tips + 1) / 2, n_tips / 2);
  
  /*Rcout << " Maximum pair score on " << n_tips << " tips: " << max_score
        << ": lg2_unrooted[n] = " << lg2_unrooted[n_tips] << " - ltms("
        << ((n_tips + 1) / 2) << ", " << (n_tips / 2) << ") = "
        << lg2_trees_matching_split((n_tips + 1) / 2, n_tips / 2) << "\n\n";*/
  
  int** score = new int*[max_splits];
  for (int i = 0; i < max_splits; i++) score[i] = new int[max_splits];
  
  uint32_t different[a.n_bins];
  int n_different, n_same, n_a_only, n_a_and_b;
  double score1, score2;
  for (int ai = 0; ai < a.n_splits; ai++) {
    for (int bi = 0; bi < b.n_splits; bi++) {
      n_different = 0;
      n_a_only = 0;
      n_a_and_b = 0;
      for (int bin = 0; bin < a.n_bins; bin++) {
        different[bin] = a.state[ai][bin] ^ b.state[bi][bin];
        n_different += count_bits_32(different[bin]);
        n_a_only += count_bits_32(a.state[ai][bin] & different[bin]);
        n_a_and_b += count_bits_32(a.state[ai][bin] & ~different[bin]);
        /*n_a_and_b += n_different - count_bits_32(a.state[ai][bin]); */
      }
      n_same = n_tips - n_different;
      /*Rcout << "  a: " << ai << ", b: " << bi << "; same = " << n_same
            << ", diff = " << n_different << "; n(a&b) = " << n_a_and_b 
            << ", aOnly = " << n_a_only << "\n";*/
      
      score1 = lg2_unrooted[n_same] - 
        lg2_trees_matching_split(n_a_and_b, n_same - n_a_and_b);
      
      score2 = lg2_unrooted[n_different] - 
        lg2_trees_matching_split(n_a_only, n_different - n_a_only);
      
      score[ai][bi] = BIG * 
        (1 - ((score1 > score2) ? score1 : score2) / max_score);
      
      /*Rcout << "    Score: 1=" << score1 << ", 2=" << score2 << ", max = "
            << ((score1 > score2) ? score1 : score2) << " = " 
            << (((score1 > score2) ? score1 : score2) / max_score) << "\n\n";*/
    }
    for (int bi = b.n_splits; bi < max_splits; bi++) {
      score[ai][bi] = BIG;
    }
  }
  for (int ai = a.n_splits; ai < max_splits; ai++) {
    for (int bi = 0; bi < max_splits; bi++) {
      score[ai][bi] = BIG;
    }
  }
  
  lap_col *rowsol = new lap_col[max_splits];
  lap_row *colsol = new lap_row[max_splits];
  cost *u = new cost[max_splits], *v = new cost[max_splits];
  
  NumericVector final_score = NumericVector::create(
    (double)((BIG * max_splits) - lap(max_splits, score, rowsol, colsol, u, v))
    * max_score / BIGL),
    final_matching (max_splits);
  
  for (int i = 0; i < max_splits; i++) {
    final_matching[i] = rowsol[i] + 1;
  }
  
  List ret = List::create(Named("score") = final_score,
                          _["matching"] = final_matching);
  
  return (ret);
}

double p_lg2_p_frac (double p) {
  return -p * log2(p);
}

double p_lg2_p (double p) {
  if (p == 0) return 0;
  if (p == 1) return 0;
  return p_lg2_p_frac(p);
}

double entropy2 (double p) {
  if (p == 0) return 0;
  if (p == 1) return 0;
  return p_lg2_p_frac(p) + p_lg2_p_frac(1 - p);
}

double entropy4 (double p1, double p2, double p3, double p4) {
  return p_lg2_p(p1) +  p_lg2_p(p2) +  p_lg2_p(p3) +  p_lg2_p(p4);
}

// [[Rcpp::export]]
List cpp_mutual_clustering (NumericMatrix x, NumericMatrix y,
                           NumericVector nTip) {
  if (x.cols() != y.cols()) {
    throw std::invalid_argument("Input splits must address same number of tips.");
  }
  SplitList a(x), b(y);
  const int max_splits = (a.n_splits > b.n_splits) ? a.n_splits : b.n_splits,
    last_bin = a.n_bins - 1,
    n_tips = nTip[0],
    unset_tips = (n_tips % 32) ? 32 - n_tips % 32 : 32;
  const uint32_t unset_mask = ~0U >> unset_tips;
  
  uint32_t b_compl[b.n_splits][b.n_bins];
  for (int i = 0; i < b.n_splits; i++) {
    for (int bin = 0; bin < last_bin; bin++) {
      b_compl[i][bin] = ~b.state[i][bin];
    }
    b_compl[i][last_bin] = b.state[i][last_bin] ^ unset_mask;
  }
  
  
  /*Rcout << " Maximum pair score on " << n_tips << " tips: " << max_score
        << ": lg2_unrooted[n] = " << lg2_unrooted[n_tips] << " - ltms("
        << ((n_tips + 1) / 2) << ", " << (n_tips / 2) << ") = "
        << lg2_trees_matching_split((n_tips + 1) / 2, n_tips / 2) << "\n\n";*/
  
  int** score = new int*[max_splits];
  for (int i = 0; i < max_splits; i++) score[i] = new int[max_splits];
  
  double a_and_b, a_and_B, A_and_b, A_and_B, 
    p1, p2;
  for (int ai = 0; ai < a.n_splits; ai++) {
    for (int bi = 0; bi < b.n_splits; bi++) {
      a_and_b = 0;
      a_and_B = 0;
      A_and_b = n_tips;
      for (int bin = 0; bin < a.n_bins; bin++) {
        a_and_b += count_bits_32(a.state[ai][bin] & b.state[bi][bin]);
        a_and_B += count_bits_32(a.state[ai][bin] & b_compl[bi][bin]);
        A_and_b -= count_bits_32(a.state[ai][bin] | b_compl[bi][bin]);
      }
      
      /* Convert to probabilities */
      a_and_b /= (double) n_tips;
      a_and_B /= (double) n_tips;
      A_and_b /= (double) n_tips;
      A_and_B = 1 - (a_and_b + a_and_B + A_and_b);
      
      p1 = a_and_b + A_and_b;
      p2 = a_and_b + a_and_B;
      
      score[ai][bi] = BIG * (1 - ((entropy2(p1) + entropy2(p2) - 
        entropy4(a_and_b, a_and_B, A_and_b, A_and_B))));
    }
    for (int bi = b.n_splits; bi < max_splits; bi++) {
      score[ai][bi] = BIG;
    }
  }
  for (int ai = a.n_splits; ai < max_splits; ai++) {
    for (int bi = 0; bi < max_splits; bi++) {
      score[ai][bi] = BIG;
    }
  }
  
  lap_col *rowsol = new lap_col[max_splits];
  lap_row *colsol = new lap_row[max_splits];
  cost *u = new cost[max_splits], *v = new cost[max_splits];
  
  NumericVector final_score = NumericVector::create(
    (double)((BIG * max_splits) - lap(max_splits, score, rowsol, colsol, u, v))
    * n_tips / BIGL),
    final_matching (max_splits);
  
  for (int i = 0; i < max_splits; i++) {
    final_matching[i] = rowsol[i] + 1;
  }
  
  List ret = List::create(Named("score") = final_score,
                          _["matching"] = final_matching);
  
  return (ret);
}

double one_overlap (int a, int b, int n) {
  if (a == b) return lg2_rooted[a] + lg2_rooted[n - a];
  if (a < b) return lg2_rooted[b] + lg2_rooted[n - a] - lg2_rooted[b - a + 1];
  return lg2_rooted[a] + lg2_rooted[n - b] - lg2_rooted[a - b + 1];
}

double one_overlap_notb (int a, int n_minus_b, int n) {
  const int b = n - n_minus_b;
  if (a == b) return lg2_rooted[b] + lg2_rooted[n_minus_b];
  if (a < b) return lg2_rooted[b] + lg2_rooted[n - a] - lg2_rooted[b - a + 1];
  return lg2_rooted[a] + lg2_rooted[n_minus_b] - lg2_rooted[a - b + 1];
}

double mpi (uint32_t* a_state, uint32_t* b_state, int n_tips, 
            int in_a, int in_b, double lg2_unrooted_n, int n_bins) {
  bool flag = true;
  
  for (int bin = 0; bin < n_bins; bin++) {
    if (a_state[bin] & b_state[bin]) {
      flag = false;
      break;
    }
  }
  if (flag) return lg2_unrooted_n - one_overlap_notb(in_a, in_b, n_tips);
  
  for (int bin = 0; bin < n_bins; bin++) {
    if ((~a_state[bin] & b_state[bin])) {
      flag = true;
      break;
    }
  }
  if (!flag) return lg2_unrooted_n - one_overlap(in_a, in_b, n_tips);
    
  for (int bin = 0; bin < n_bins; bin++) {
    if ((a_state[bin] & ~b_state[bin])) {
      flag = false;
      break;
    }
  }
  if (flag) return lg2_unrooted_n - one_overlap(in_a, in_b, n_tips);
  
  for (int bin = 0; bin < n_bins; bin++) {
    if (~(a_state[bin] | b_state[bin])) {
      flag = true;
      break;
    }
  }
  if (!flag) return lg2_unrooted_n - one_overlap_notb(in_a, in_b, n_tips);
  
  return 0;
}

// [[Rcpp::export]]
List cpp_mutual_phylo (NumericMatrix x, NumericMatrix y,
                              NumericVector nTip) {
  if (x.cols() != y.cols()) {
    throw std::invalid_argument("Input splits must address same number of tips.");
  }
  SplitList a(x), b(y);
  const int max_splits = (a.n_splits > b.n_splits) ? a.n_splits : b.n_splits,
    n_tips = nTip[0];
  const double lg2_unrooted_n = lg2_unrooted[n_tips],
    max_score = lg2_unrooted_n - one_overlap((n_tips + 1) / 2, n_tips / 2, n_tips);
  int in_a[a.n_splits], in_b[b.n_splits];
  
  for (int i = 0; i < a.n_splits; i++) {
    in_a[i] = 0;
    for (int bin = 0; bin < a.n_bins; bin++) {
      in_a[i] += count_bits_32(a.state[i][bin]);
    }
  }
  for (int i = 0; i < b.n_splits; i++) {
    in_b[i] = 0;
    for (int bin = 0; bin < b.n_bins; bin++) {
      in_b[i] += count_bits_32(b.state[i][bin]);
    }
  }
  
  int** score = new int*[max_splits];
  for (int i = 0; i < max_splits; i++) score[i] = new int[max_splits];
  
  for (int ai = 0; ai < a.n_splits; ai++) {
    for (int bi = 0; bi < b.n_splits; bi++) {
      score[ai][bi] = BIG * (1 - 
        (mpi(a.state[ai], b.state[bi], n_tips, in_a[ai], in_b[bi],
            lg2_unrooted_n, a.n_bins) / max_score));
    }
    for (int bi = b.n_splits; bi < max_splits; bi++) {
      score[ai][bi] = BIG;
    }
  }
  for (int ai = a.n_splits; ai < max_splits; ai++) {
    for (int bi = 0; bi < max_splits; bi++) {
      score[ai][bi] = BIG;
    }
  }
  
  lap_col *rowsol = new lap_col[max_splits];
  lap_row *colsol = new lap_row[max_splits];
  cost *u = new cost[max_splits], *v = new cost[max_splits];
  
  NumericVector final_score = NumericVector::create(
    (double)((BIG * max_splits) - lap(max_splits, score, rowsol, colsol, u, v))
    * max_score / BIGL),
    final_matching (max_splits);
  
  for (int i = 0; i < max_splits; i++) {
    final_matching[i] = rowsol[i] + 1;
  }
  
  List ret = List::create(Named("score") = final_score,
                          _["matching"] = final_matching);
  
  return (ret);
}
