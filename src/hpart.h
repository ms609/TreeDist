// src/hpart.h
#pragma once
#include <vector>
#include <cstdint>
#include <TreeTools/assert.h> // for ASSERT
#include <Rcpp.h>

namespace TreeDist {

// Lookup table size - adjust based on your typical x values
constexpr size_t LOOKUP_SIZE = 8192;

// Pre-computed lookup table
constexpr auto make_x_log_x_table() {
  std::array<double, LOOKUP_SIZE> table{};
  table[0] = 0.0;
  table[1] = 0.0;
  for (size_t i = 2; i < LOOKUP_SIZE; ++i) {
    table[i] = i * std::log(static_cast<double>(i));
  }
  return table;
}

// Static lookup table initialized at compile time
static constexpr auto x_log_x_table = make_x_log_x_table();

static inline double x_log_x(size_t x) {
  if (x < LOOKUP_SIZE) [[likely]] {
    return x_log_x_table[x];
  } else {
    // Fallback for large x values
    return x * std::log(static_cast<double>(x));
  }
}
// Alternative version with branch-free lookup for even better performance
// when you're certain x will usually be small
static inline double x_log_x_branchless(size_t x) {
  // Use lookup if x < LOOKUP_SIZE, otherwise compute
  const bool use_table = x < LOOKUP_SIZE;
  const double table_result = x < LOOKUP_SIZE ? x_log_x_table[x] : 0.0;
  const double computed_result = x * std::log(static_cast<double>(x));
  
  // Branchless selection (compiler will optimize this well)
  return use_table ? table_result : computed_result;
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
