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

struct HNode {
  std::vector<size_t> children;   // indices of children in HPart.nodes
  int label = -1;                 // for tips; counting from zero
  uint64_t* bitset;               // pointer into leaf set pool
                                  // Faster on the heap than the stack.
  
  int leaf_count = 0;
  double x_log_x = 0;             // where x = leaf_count
  bool all_kids_leaves = true;
  int n_tip;
  size_t n_block;
  double n_tip_reciprocal;
};

struct HPart {
  std::vector<HNode> nodes;
  std::vector<uint64_t> bitset_pool;
  size_t n_tip;
  size_t n_block;
  
  size_t root;
  
  double entropy = std::numeric_limits<double>::min();
  double hier_entropy = std::numeric_limits<double>::min();
  
  // Constructor for building from n_tip and number of nodes
  HPart(size_t n_tip_, size_t vec_size)
    : nodes(vec_size), n_tip(n_tip_)
  {
    n_block = (n_tip + 63) / 64;
    const double n_tip_recip = 1 / static_cast<double>(n_tip);
    
    bitset_pool.assign((vec_size - 1) * n_block, 0ULL);
    
    for (size_t i = 1; i < vec_size; ++i) {
      TreeDist::HNode& n = nodes[i];
      n.n_tip = n_tip;
      n.n_block = n_block;
      n.n_tip_reciprocal = n_tip_recip;
      n.bitset = &bitset_pool[(i - 1) * n_block];
    }
  }
  
  HPart(const HPart& other) {
    // Copy pool and metadata
    bitset_pool = other.bitset_pool;
    nodes = other.nodes;
    n_tip = other.n_tip;
    n_block = other.n_block;
    
    root = other.root;
    
    entropy = other.entropy;
    hier_entropy = other.hier_entropy;
    
    // Rewire pointers
    for (size_t i = 0; i < nodes.size(); ++i) {
      nodes[i].bitset = bitset_pool.data() + (i * n_block);
    }
  }
  
  HPart& operator=(const HPart& other) {
    if (this != &other) {
      nodes = other.nodes;
      bitset_pool = other.bitset_pool;
      n_tip = other.n_tip;
      n_block = other.n_block;
      
      root = other.root;
      
      entropy = other.entropy;
      hier_entropy = other.hier_entropy;
      
      for (size_t i = 0; i < nodes.size(); ++i) {
        nodes[i].bitset = bitset_pool.data() + (i * n_block);
      }
    }
    return *this;
  }
};
}

SEXP clone_hpart(SEXP hpart_ptr);
void relabel_hpart(SEXP hpart_ptr, const std::vector<int>& map);
