#include <Rcpp.h>
using namespace Rcpp;
#include <stdint.h>

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
  static int n_tips, n_splits, n_bins;
  uint32_t state[3200][100]; /* Maximum tips supported: 3200 */
  SplitList(LogicalMatrix);
  int n() {return n_splits;}
  int bins() {return n_bins;}
};

SplitList::SplitList(LogicalMatrix x) {
  n_tips = x.rows();
  n_splits = x.cols();
  if (n_tips < 1) throw std::range_error("No tips present.");
  if (n_splits < 1) throw std::range_error("No splits present.");
  if (n_tips > 3200) {
    throw std::range_error("No more than 3200 tips can be supported. Please contact the maintainer if you need to use more!");
  }
  n_bins = (n_tips - 1) / 32 + 1;
  int current_block;
  
  for (int i = 0; i < n_splits; i++) {
    current_block = -1;
    for (int j = 0; j < n_tips; j++) {
      if (j % 32 == 0) ++current_block;
      state[i][current_block] <<= 1;
      if (x(j, i)) state[i][current_block]++;
    }
  }
}

// [[Rcpp::export]]
IntegerMatrix matching_split_distance (LogicalMatrix x, LogicalMatrix y) {
  if (x.rows() != y.rows()) {
    throw std::range_error("Input matrices must contain same number of rows.");
  }
  SplitList splits_x = SplitList(x);
  SplitList splits_y = SplitList(y);
  IntegerMatrix score(splits_x.n_splits, splits_y.n_splits);
  static int n_tips = splits_x.n_tips,
    half_tips = n_tips / 2;
  
  for (int xi = 0; xi < x.cols(); xi++) {
    for (int yi = 0; yi < y.cols(); yi++) {
      score(xi, yi) = 0;
      for (int bin = 0; bin < splits_x.bins(); bin++) {
        score(xi, yi) += count_bits_32(splits_x.state[xi][bin] ^ splits_y.state[yi][bin]);
      }
      if (score(xi, yi) < half_tips) score(xi, yi) = n_tips - score(xi, yi);
    }
  }
  
  return (score);
}
