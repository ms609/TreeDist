// src/hpart.h
#pragma once
#include <vector>
#include <cstdint>
#include <TreeTools/assert.h> // for ASSERT
#include <Rcpp.h>

namespace TreeDist {
struct HNode {
  std::vector<size_t> children;   // indices of children in HPart.nodes
  int label = -1;                 // for tips; counting from zero
  std::vector<uint64_t> bitset;   // leaf set
  int leaf_count = 0;
  bool all_kids_leaves = true;
  int n_tip = 0;
};

struct HPart {
  std::vector<HNode> nodes;  // owns all nodes
  size_t root;
  double entropy = std::numeric_limits<double>::min();
};
}

SEXP clone_hpart(SEXP hpart_ptr);
void relabel_hpart(SEXP hpart_ptr, const std::vector<int>& map);
