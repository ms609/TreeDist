#include <TreeTools/SplitList.h>
#include <Rcpp/Lightest>

// Provide the MCI table definitions and implementation in this TU.
#define TREEDIST_MCI_IMPLEMENTATION
#include <TreeDist/mutual_clustering_impl.h>

#include "tree_distances.h"

// Populate lookup tables at library load time.
__attribute__((constructor))
  void initialize_ldf() {
    TreeDist::init_lg2_tables(SL_MAX_TIPS);
  }

// Thin wrapper that exercises the installable-header version of
// mutual_clustering_score(), for test coverage and downstream validation.
// [[Rcpp::export]]
double cpp_mci_impl_score(const Rcpp::RawMatrix& x,
                          const Rcpp::RawMatrix& y,
                          int n_tips) {
  using TreeTools::SplitList;

  const SplitList a(x);
  const SplitList b(y);
  TreeDist::LapScratch scratch;

  // Build arrays matching the header's raw-pointer API types.
  std::vector<const splitbit*> a_ptrs(a.n_splits);
  std::vector<const splitbit*> b_ptrs(b.n_splits);
  std::vector<TreeDist::split_int> a_in(a.n_splits);
  std::vector<TreeDist::split_int> b_in(b.n_splits);
  for (TreeDist::int32 i = 0; i < a.n_splits; ++i) {
    a_ptrs[i] = a.state[i];
    a_in[i] = static_cast<TreeDist::split_int>(a.in_split[i]);
  }
  for (TreeDist::int32 i = 0; i < b.n_splits; ++i) {
    b_ptrs[i] = b.state[i];
    b_in[i] = static_cast<TreeDist::split_int>(b.in_split[i]);
  }

  return TreeDist::mutual_clustering_score(
    a_ptrs.data(), a_in.data(), a.n_splits,
    b_ptrs.data(), b_in.data(), b.n_splits,
    a.n_bins, static_cast<TreeDist::int32>(n_tips),
    scratch);
}

// [[Rcpp::export]]
int cpp_max_tips() {
  return static_cast<int>(SL_MAX_TIPS);
}
