#include "hpart.h"
#include <Rcpp.h>

// [[Rcpp::depends(TreeTools)]]
#include <TreeTools/renumber_tree.h> // for preorder_edges_and_nodes

using namespace Rcpp;
using Node = TreeDist::HNode;

// [[Rcpp::export]]
SEXP build_hpart_from_phylo(List phy) {
  IntegerMatrix edge = phy["edge"];
  CharacterVector tip_label = phy["tip.label"];
  int n_tip = tip_label.size();
  int n_node = phy["Nnode"];
  
  IntegerMatrix reordered = TreeTools::preorder_edges_and_nodes(edge(_,0), edge(_,1));
  
  const size_t vec_size = n_tip + n_node + 1;
  std::vector<std::vector<size_t>> children(vec_size);
  for (size_t i = 0; i < children.size(); ++i) {
    children[i].reserve(2);
  }
  for (size_t i = 0; i < reordered.nrow(); ++i) {
    size_t p = reordered(i, 0);
    size_t c = reordered(i, 1);
    children[p].push_back(c);
  }
  
  HPart* hpart = new HPart();
  hpart->nodes.resize(vec_size);
  
  // Initialize all nodes to empty
  const int n_block = (n_tip + 63) / 64;
  for (size_t i = 1; i < vec_size; ++i) {
    hpart->nodes[i].nTip = nTip;
    hpart->nodes[i].bitset.resize(nBlock, 0);
  }
  
  // Initialize tips
  for (size_t i = 1; i <= n_tip; ++i) {
    const size_t bit_index = i - 1; // 0-based indexing
    const size_t vector_pos = bit_index / 64; 
    const size_t bit_pos_in_block = bit_index % 64;
    hpart->nodes[i].bitset[vector_pos] = 1ULL << bit_pos_in_block;
    hpart->nodes[i].leafCount = 1;
    hpart->nodes[i].calc_entropy();
  }
  
  // Traverse nodes in postorder
  for (size_t i = vec_size; i > (size_t)n_tip; --i) {
    hpart->nodes[i].children.reserve(children[i].size());
    
    for (size_t child_id : children[i]) {
      auto& child_node = &hpart->nodes[child_id];
      
      hpart->nodes[i].children.push_back(child_node);
      const size_t child_leaves = child_node.leaf_count;
      if (child_leaves > 1) {
        hpart->nodes[i].all_kids_leaves = false;
      }
      hpart->nodes[i].leaf_count += child_leaves;
      
      for (size_t chunk = 0; chunk < nodes[i].bitset.size(); ++chunk) {
        hpart->nodes[i].bitset[chunk] |= child_node.bitset[chunk];
      }
    }
    hpart->nodes[i].calc_entropy();
  }
  
  hpart->root = hpart->nodes[n_tip + 1];
  
  return Rcpp::XPtr<HPart>(hpart, true);
}
