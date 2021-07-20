#ifndef _TREEDIST_SPLITROSTER_HPP
#define _TREEDIST_SPLITROSTER_HPP

#include <stdint.h>
#include <Rcpp.h>
#include "tree_distances.h"
#include "SplitList.h"

using namespace Rcpp;

class SplitRoster {
public:
  int16 n_splits, n_trees;
};

#endif
