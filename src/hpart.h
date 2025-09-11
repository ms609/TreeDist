// src/hpart.h
#pragma once
#include <vector>
#include <cstdint>
#include <TreeTools/assert.h> // for ASSERT
#include <Rcpp.h>

namespace TreeDist {

static inline double x_log_x(size_t x) {
  return x > 1 ? x * std::log(x) : 0.0;
}

struct HNode {
  std::vector<size_t> children;   // indices of children in HPart.nodes
  int label = -1;                 // for tips; counting from zero
  std::vector<uint64_t> bitset;   // leaf set
  int leaf_count = 0;
  double x_log_x = 0;             // where x = leaf_count
  bool all_kids_leaves = true;
  int n_tip;
  double n_tip_reciprocal;
};

struct HPart {
  std::vector<HNode> nodes;
  size_t root;
  double entropy = std::numeric_limits<double>::min();
  double hier_entropy = std::numeric_limits<double>::min();
};
}

SEXP clone_hpart(SEXP hpart_ptr);
void relabel_hpart(SEXP hpart_ptr, const std::vector<int>& map);
