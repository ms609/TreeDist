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
  int n_tips, n_splits, n_bins;
  uint32_t state[3200][100]; /* Maximum tips supported: 3200 */
  SplitList(LogicalMatrix);
  int n() {return n_splits;}
  int bins() {return n_bins;}
  int tips() {return n_tips;}
};

SplitList::SplitList(LogicalMatrix x) {
  int current_block;
  
  n_tips = x.rows();
  n_splits = x.cols();
  
  if (n_tips < 1) throw std::invalid_argument("No tips present.");
  if (n_splits < 1) throw std::invalid_argument("No splits present.");
  if (n_tips > 3200) {
    throw std::length_error("No more than 3200 tips can be supported. Please contact the TreeDist maintainer if you need to use more!");
  }
  
  n_bins = (n_tips - 1) / 32 + 1;
  
  for (int i = 0; i < n_splits; i++) {
    current_block = -1;
    for (int j = 0; j < n_tips; j++) {
      if (j % 32 == 0) {
        ++current_block;
        state[i][current_block] = 0;
      }
      state[i][current_block] <<= 1;
      if (x(j, i)) state[i][current_block]++;
    }
  }
}

// [[Rcpp::export]]
List cpp_matching_split_distance (LogicalMatrix x, LogicalMatrix y) {
  if (x.rows() != y.rows()) {
    throw std::invalid_argument("Input matrices must contain same number of rows.");
  }
  SplitList a(x);
  SplitList b(y);
  const int max_splits = (a.n_splits > b.n_splits) ? a.n_splits : b.n_splits;
  const int split_diff = max_splits -
    ((a.n_splits > b.n_splits) ? b.n_splits : a.n_splits);
  int** score = new int*[max_splits];
  for (int i = 0; i < max_splits; i++) score[i] = new int[max_splits];
  const int n_tips = a.n_tips, half_tips = n_tips / 2;
  
  /*Rcout << "Working over " << a.n() << " (" << a.n_splits << ", " << x.cols() 
        << ") and " << b.n() << " (" << b.n_splits << ", " << y.cols() 
        << ") splits.\n\n";*/
  
  for (int ai = 0; ai < a.n_splits; ai++) {
    for (int bi = 0; bi < b.n_splits; bi++) {
      score[ai][bi] = 0;
      for (int bin = 0; bin < a.bins(); bin++) {
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

