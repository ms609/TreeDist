// src/hpart.h
#pragma once
#include <vector>
#include <cstdint>
#include <Rcpp.h>

struct Node {
  std::vector<int> children;
  int labelIndex = -1;                  // for tips
  std::vector<uint64_t> bitset;         // leaf set
  int leafCount = 0;
  bool allKidsLeaves = true;
};

struct HPart {
  std::vector<Node> nodes;
  int nTips = 0;
  int nBlocks = 0;
  int rootIndex = -1;
};
