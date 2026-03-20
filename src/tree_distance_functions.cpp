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
