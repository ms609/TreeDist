#include "hpart.h"
#include <Rcpp.h>

// [[Rcpp::depends(TreeTools)]]
#include <TreeTools/renumber_tree.h> // for preorder_edges_and_nodes

using namespace Rcpp;

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
  for (size_t i = 0; i < (size_t)reordered.nrow(); ++i) {
    size_t p = reordered(i, 0);
    size_t c = reordered(i, 1);
    children[p].push_back(c);
  }
  
  TreeDist::HPart* hpart = new TreeDist::HPart();
  hpart->nodes.resize(vec_size);
  
  // Initialize all nodes to empty
  const int n_block = (n_tip + 63) / 64;
  for (size_t i = 1; i < vec_size; ++i) {
    hpart->nodes[i].n_tip = n_tip;
    hpart->nodes[i].bitset.resize(n_block, 0);
  }
  
  // Initialize tips
  for (size_t i = 1; i <= (size_t)n_tip; ++i) {
    const size_t bit_index = i - 1; // 0-based indexing
    const size_t vector_pos = bit_index / 64; 
    const size_t bit_pos_in_block = bit_index % 64;
    auto &node_i = hpart->nodes[i];
    node_i.bitset[vector_pos] = 1ULL << bit_pos_in_block;
    node_i.leaf_count = 1;
    node_i.label = i - 1;
    node_i.calc_entropy();
  }
  
  // Traverse nodes in postorder
  for (size_t i = vec_size - 1; i > (size_t)n_tip; --i) {
    auto &node_i = hpart->nodes[i];
    node_i.children.reserve(children[i].size());
    
    for (size_t child_id : children[i]) {
      const auto child_node = &hpart->nodes[child_id];
      
      node_i.children.push_back(child_id);
      const size_t child_leaves = child_node->leaf_count;
      if (child_leaves > 1) {
        node_i.all_kids_leaves = false;
      }
      node_i.leaf_count += child_leaves;
      
      for (size_t chunk = 0; chunk < node_i.bitset.size(); ++chunk) {
        node_i.bitset[chunk] |= child_node->bitset[chunk];
      }
    }
    node_i.calc_entropy();
  }
  
  hpart->root = n_tip + 1;
  
  return Rcpp::XPtr<TreeDist::HPart>(hpart, true);
}

// helper: get index of a node pointer within hpart->nodes
inline size_t node_index(const TreeDist::HNode* node,
                         const std::vector<TreeDist::HNode>& nodes) {
  return static_cast<size_t>(node - &nodes[0]);
}

// [[Rcpp::export]]
Rcpp::IntegerMatrix hpart_to_edge(SEXP hpart_xptr) {
  Rcpp::XPtr<TreeDist::HPart> hpart(hpart_xptr);
  TreeDist::HPart* hp = hpart.get();
  
  int n_tip = hp->nodes[1].n_tip;
  int n_node = static_cast<int>(hp->nodes.size()) - n_tip - 1;
  int n_edge = n_tip + n_node - 1;
  
  // Assign IDs: tips get their label, internal nodes get sequential IDs
  std::vector<int> node_ids(hp->nodes.size(), -1);
  int next_index = n_tip + 1;
  for (int i = 1; i <= n_tip; ++i) {
    node_ids[i] = hp->nodes[i].label + 1; // 1-based R tip index
  }
  for (size_t i = n_tip + 1; i < hp->nodes.size(); ++i) {
    node_ids[i] = next_index++;
  }
  
  // Collect edges
  std::vector<std::pair<int,int>> edges;
  edges.reserve(n_edge);
  
  for (size_t i = n_tip + 1; i < hp->nodes.size(); ++i) {
    int parent_id = node_ids[i];
    for (size_t cidx : hp->nodes[i].children) {
      int child_id = node_ids[cidx];
      edges.emplace_back(parent_id, child_id);
    }
  }
  
  // Build R matrix
  IntegerMatrix edge_mat(edges.size(), 2);
  for (size_t i = 0; i < edges.size(); ++i) {
    edge_mat(i, 0) = edges[i].first;
    edge_mat(i, 1) = edges[i].second;
  }
  
  return edge_mat;
}

// [[Rcpp::export]]
SEXP clone_hpart(SEXP hpart_ptr) {
  Rcpp::XPtr<TreeDist::HPart> src(hpart_ptr);
  
  TreeDist::HPart* copy = new TreeDist::HPart(*src);
  
  return Rcpp::XPtr<TreeDist::HPart>(copy, true);
}
