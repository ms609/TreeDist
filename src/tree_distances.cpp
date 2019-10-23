#include <Rcpp.h>
using namespace Rcpp;
#include <stdint.h>
#include "SplitList.h"
#include "lap.h"

uint32_t bitcounts[65536]; // the bytes representing bit count of each number 0-65535
__attribute__((constructor))
  void initialize_bitcounts()
  {
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


// [[Rcpp::export]]
List cpp_matching_split_distance (NumericMatrix x, NumericMatrix y, 
                                  NumericVector nTip) {
  if (x.cols() != y.cols()) {
    throw std::invalid_argument("Input splits must address same number of tips.");
  }
  SplitList a(x);
  SplitList b(y);
  const int max_splits = (a.n_splits > b.n_splits) ? a.n_splits : b.n_splits;
  const int split_diff = max_splits -
    ((a.n_splits > b.n_splits) ? b.n_splits : a.n_splits);
  int** score = new int*[max_splits];
  for (int i = 0; i < max_splits; i++) score[i] = new int[max_splits];
  const int n_tips = nTip[0], half_tips = n_tips / 2;
  
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
  cost *u = new cost[max_splits], *v = new cost[max_splits];
  
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
List cpp_nye_distance (NumericMatrix x, NumericMatrix y,
                       NumericVector nTip) {
  if (x.cols() != y.cols()) {
    throw std::invalid_argument("Input splits must address same number of tips.");
  }
  SplitList a(x);
  SplitList b(y);
  const int max_splits = (a.n_splits > b.n_splits) ? a.n_splits : b.n_splits;
  
  int a_and_b, a_and_B, A_and_b, A_and_B,
      a_or_b,  a_or_B,  A_or_b,  A_or_B;
  double ars_ab, ars_aB, ars_Ab, ars_AB,
         min_ars_both, min_ars_either;
  
  
  int** score = new int*[max_splits];
  for (int i = 0; i < max_splits; i++) score[i] = new int[max_splits];
  const int n_tips = nTip[0];
  
  Rcout << "Working over " << a.n_splits << " (" << a.n_splits << ", " << x.rows() 
          << ") and " << b.n_splits << " (" << b.n_splits << ", " << y.rows() 
          << ") splits.\n\n";
  
  for (int ai = 0; ai < a.n_splits; ai++) {
    for (int bi = 0; bi < b.n_splits; bi++) {
      a_and_b = 0;
      a_and_B = 0;
      A_and_b = 0;
      a_or_b = 0; /* na_and_nb = n - a_or_b */
      for (int bin = 0; bin < a.n_bins; bin++) {
        a_and_b += count_bits_32( a.state[ai][bin] &  b.state[bi][bin]);
        a_and_B += count_bits_32( a.state[ai][bin] & ~b.state[bi][bin]);
        A_and_b += count_bits_32(~a.state[ai][bin] &  b.state[bi][bin]);
        a_or_b  += count_bits_32( a.state[ai][bin] |  b.state[bi][bin]);
        
        Rcout << "\n- x = " << ai << ", y = " << bi << ", bin " << bin << ": A = "
              << a.state[ai][bin] << ", B = " << b.state[bi][bin] << ".\nA & B = "
              << count_bits_32( a.state[ai][bin] &  b.state[bi][bin]) 
              << "; A & ~B = " << count_bits_32( a.state[ai][bin] & ~b.state[bi][bin])
              << " (" << a_and_B << ") "
              << "; ~A & B = " << count_bits_32(~a.state[ai][bin] &  b.state[bi][bin])
              << "; A | B = " << count_bits_32( a.state[ai][bin] |  b.state[bi][bin])
              <<  " (" << n_tips << " tips).\n";
      }
      A_and_B = n_tips - a_or_b;
      A_or_B  = n_tips - a_and_b;
      a_or_B  = n_tips - A_and_b;
      A_or_b  = n_tips - a_and_B;
      
      ars_ab = (double) a_and_b / (double) a_or_b;
      ars_Ab = (double) A_and_b / (double) A_or_b;
      ars_aB = (double) a_and_B / (double) a_or_B;
      ars_AB = (double) A_and_B / (double) A_or_B;
      
      min_ars_both = (ars_ab < ars_AB) ? ars_ab : ars_AB;
      min_ars_either = (ars_aB < ars_Ab) ? ars_aB : ars_Ab;
      
      /* LAP will look to minimize an integer. max(ars) is between 0 and 1. */
      score[ai][bi] = BIG - (BIG * ((min_ars_both > min_ars_either) ? 
                                     min_ars_both : min_ars_either));
      Rcout /*<< "  ab = " << ars_ab << "; AB = " <<  ars_AB
            << "\n  Ab = " << ars_Ab << " aB = " << ars_aB
            << "\n  Min Ars both = " << min_ars_both 
            << " Min Ars Either = " << min_ars_either */
            << " Score: " << score[ai][bi] << "\n";
      
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
  
  Rcout << "\n * LAP solution = " 
        << lap(max_splits, score, rowsol, colsol, u, v) << "\n";
  Rcout << "(" << (BIG * max_splits) << " - " 
        << lap(max_splits, score, rowsol, colsol, u, v) << ") / " 
        << BIG << "\n";
  NumericVector final_score = NumericVector::create(
    (double)((BIG * max_splits) - lap(max_splits, score, rowsol, colsol, u, v))
    / (double) BIG),
    final_matching (max_splits);
  
  for (int i = 0; i < max_splits; i++) {
    final_matching[i] = rowsol[i] + 1;
  }
  
  List ret = List::create(Named("score") = final_score,
                          _["matching"] = final_matching);
  
  return (ret);
}
