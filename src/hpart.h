// src/hpart.h
#pragma once
#include <vector>
#include <cstdint>
#include <Rcpp.h>

namespace TreeDist {
struct HNode {
  std::vector<HNode*> children;
  int labelIndex = -1;                  // for tips
  std::vector<uint64_t> bitset;         // leaf set
  int leafCount = 0;
  bool allKidsLeaves = true;
  int nTip = 0;
};

struct HPart {
  std::vector<HNode> nodes;  // owns all nodes
  HNode* root = nullptr;     // pointer into nodes
};
}
