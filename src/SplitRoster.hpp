#ifndef _TREEDIST_SPLITROSTER_HPP
#define _TREEDIST_SPLITROSTER_HPP

#include <stdint.h>
#include <Rcpp.h>
#include "ints.hpp"
#include "tree_distances.hpp"
#include "SplitList.hpp"

using namespace Rcpp;

class SplitRoster {
  std::unique_ptr<SplitList[]> splits;
  std::unique_ptr<SplitList[]> roster;
  int16 n_bins;
public:
  int16 n_tips, n_splits, n_trees;
};

#endif
