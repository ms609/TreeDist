// src/hpart.h
#pragma once
#include <assert>
#include <vector>
#include <cstdint>
#include <Rcpp.h>

namespace TreeDist {
struct HNode {
  std::vector<HNode*> children;
  int label_index = -1;                  // for tips
  std::vector<uint64_t> bitset;         // leaf set
  int leaf_count = 0;
  bool all_kids_leaves = true;
  int n_tip = 0;
  double entropy = 0;
  
  void calc_entropy() {
    std::assert(this->leaf_count > 0);
    std::assert(this->n_tip > 0);
    double p = static_cast<double>(this->leaf_count) / this->n_tip;
    this->entropy = -p * std::log(p);
  }
};

struct HPart {
  std::vector<HNode> nodes;  // owns all nodes
  HNode* root = nullptr;     // pointer into nodes
};
}
