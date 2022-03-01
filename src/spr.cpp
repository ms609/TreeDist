#include <Rcpp/Lightest>
#include <TreeTools/assert.h> // for ASSERT
#include <TreeTools/types.h> // for intx
#include <TreeTools/SplitList.h> // for SplitList, count_bits

using namespace Rcpp;

// Equivalent to TREETOOLS_SPLITLIST_INIT
// Is it possible to invoke the constructor in a more elegant way,
// and to avoid duplicating the process from TreeTools?
__attribute__((constructor))                                     \
  void _treedist_initialize_bitcounts() {                        \
    for (int i = 65536; i--; ) {                                 \
      int16 n_bits = 0;                                          \
      for (int j = 16; j--; ) {                                  \
        if (i & (1 << j)) n_bits += 1;                           \
      }                                                          \
      TreeTools::bitcounts[i] = n_bits;                          \
    }                                                            \
  }   
  
// [[Rcpp::export]]
IntegerVector mismatch_size (const RawMatrix x, const RawMatrix y) {
  // Rcout << "\n Debugging mismatch_size()\n";
  const int16 n_split = x.rows();
  if (n_split != y.rows()) {
    throw std::invalid_argument("Input splits contain same number of splits.");
  }
  if (!x.hasAttribute("nTip")) {
    Rcpp::stop("`x` lacks nTip attribute");
  }
  if (!y.hasAttribute("nTip")) {
    Rcpp::stop("`y` lacks nTip attribute");
  }
  const int16 n_tip = x.attr("nTip");
  if (n_tip != int16(y.attr("nTip"))) {
    Rcpp::stop("`x` and `y` differ in `nTip`");
  }
  
  const TreeTools::SplitList a(x), b(y);
  const int16
    last_bin = a.n_bins - 1,
    unset_tips = (n_tip % SL_BIN_SIZE) ? SL_BIN_SIZE - n_tip % SL_BIN_SIZE : 0
  ;
  const splitbit all_ones = ~(splitbit(0U));
  const splitbit unset_mask = all_ones >> unset_tips;

  int16 ret_ptr = n_split * n_split;
  IntegerVector ret(ret_ptr);
  for (int16 bi = b.n_splits; bi--; ) {
    // Rcout << "a = " << ai << ".\n";
    for (int16 ai = a.n_splits; ai--; ) {
      // Rcout << "  - b = " << bi << ".\n";
      
      --ret_ptr;
      
      // Rcout << "    - last_bin: " << ((a.state[ai][last_bin] ^ b.state[bi][last_bin]) & unset_mask)
      //       << " = " << TreeTools::count_bits(
      // (a.state[ai][last_bin] ^ b.state[bi][last_bin]) & unset_mask
      //       ) << "\n";
      ret[ret_ptr] = TreeTools::count_bits(
        (a.state[ai][last_bin] ^ b.state[bi][last_bin]) & unset_mask
        );
      for (int16 bin = last_bin; bin--; ) {
        // Rcout << "    - bin = " << bin << ".\n";
        // Rcout << "      " << (a.state[ai][bin]);
        // Rcout << " ^ " << (b.state[bi][bin]);
        // Rcout << " = " << (a.state[ai][bin] ^ b.state[bi][bin]) << std::endl;
        ret[ret_ptr] += TreeTools::count_bits(a.state[ai][bin] ^ b.state[bi][bin]);
        // Rcout << "      ret[" << ret_ptr << "] = " 
        //       << TreeTools::count_bits(a.state[ai][bin] ^ b.state[bi][bin]) 
        //       <<".\n";
      }
    }
  }
  return ret;
}


// [[Rcpp::export]]
IntegerVector confusion (const RawMatrix x, const RawMatrix y) {
  const int16 n_split = x.rows();
  if (n_split != y.rows()) {
    throw std::invalid_argument("Input splits contain same number of splits.");
  }
  if (!x.hasAttribute("nTip")) {
    Rcpp::stop("`x` lacks nTip attribute");
  }
  if (!y.hasAttribute("nTip")) {
    Rcpp::stop("`y` lacks nTip attribute");
  }
  const int16 n_tip = x.attr("nTip");
  if (n_tip != int16(y.attr("nTip"))) {
    Rcpp::stop("`x` and `y` differ in `nTip`");
  }
  
  const TreeTools::SplitList a(x), b(y);
  const int16
    n_bin = a.n_bins,
    confusion_size = 4
  ;
  IntegerVector ret(n_split * n_split * confusion_size);
  int *ret_ptr = ret.end();
  for (int16 bi = n_split; bi--; ) {
    const int16
      nb = b.in_split[bi],
      nB = n_tip - nb
    ;
    
    for (int16 ai = n_split; ai--; ) {
      
      // x divides tips into a|A; y divides tips into b|B
      int16 a_and_b = 0;
      for (int16 bin = n_bin; bin--; ) {
        a_and_b += TreeTools::count_bits(a.state[ai][bin] & b.state[bi][bin]);
      }
      
      const int16
        na = a.in_split[ai],
        a_and_B = na - a_and_b,
        A_and_b = nb - a_and_b,
        A_and_B = nB - a_and_B
      ;
      *(--ret_ptr) = A_and_B;
      *(--ret_ptr) = A_and_b;
      *(--ret_ptr) = a_and_B;
      *(--ret_ptr) = a_and_b;
    }
  }
  ret.attr("dim") = Dimension(confusion_size, n_split, n_split);
  return ret;
}

