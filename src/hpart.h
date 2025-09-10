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
  double entropy = 0;
  
  void calc_entropy() {
    ASSERT(this->leaf_count > 0);
    ASSERT(this->n_tip > 0);
    double p = static_cast<double>(this->leaf_count) / this->n_tip;
    this->entropy = -p * std::log(p);
  }
};

struct HPart {
  std::vector<HNode> nodes;  // owns all nodes
  size_t root;
};
}

SEXP clone_hpart(SEXP hpart_ptr);
void relabel_hpart(SEXP hpart_ptr, const std::vector<int>& map);
