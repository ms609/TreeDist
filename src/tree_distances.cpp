#include <Rcpp.h>
using namespace Rcpp;
#include <stdint.h>
#include "lap.h"

const uint32_t right16bits = 65535;
const uint32_t powers_of_two[16] = {1, 2, 4, 8, 16, 32, 64, 128, 256, 512,
                                    1024, 2048, 4096, 8192, 16384, 32768};

uint32_t bitcounts[65536]; // the bytes representing bit count of each number 0-65535
__attribute__((constructor))
  void initialize_array()
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

class SplitList {
public:
  int n_splits, n_bins;
  uint32_t state[3200][100]; /* Maximum tips supported: 3200 */
  SplitList(NumericMatrix);
};

SplitList::SplitList(NumericMatrix x) {
  n_splits = x.rows();
  n_bins = x.cols();
  
  if (n_bins < 1) throw std::invalid_argument("No tips present.");
  if (n_splits < 1) throw std::invalid_argument("No splits present.");
  if (n_bins > 100) {
    throw std::length_error("No more than 3200 tips can be supported. Please contact the TreeDist maintainer if you need to use more!");
  }
  
  for (int i = 0; i < n_splits; i++) {
    for (int j = 0; j < n_bins; j++) {
      state[i][j] = (uint32_t) x(i, j);
    }
  }
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
    final_matching[i] = rowsol[i];
  }
  
  List ret = List::create(Named("score") = final_score,
                          _["matching"] = final_matching);
  
  return (ret);
}

