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
  }
  
  // Traverse nodes in postorder
  for (size_t i = vec_size - 1; i > (size_t)n_tip; --i) {
    auto &node_i = hpart->nodes[i];
    node_i.children.reserve(children[i].size());
    
    for (const size_t child_id : children[i]) {
      const auto child_node = &hpart->nodes[child_id];
      
      node_i.children.push_back(child_id);
      const size_t child_leaves = child_node->leaf_count;
      if (child_leaves > 1) {
        node_i.all_kids_leaves = false;
      }
      node_i.leaf_count += child_leaves;
      
      std::transform(
        node_i.bitset.begin(), node_i.bitset.end(),
        child_node->bitset.begin(),
        node_i.bitset.begin(),
        std::bit_or<uint64_t>()
      );
    }
    node_i.x_log_x = TreeDist::x_log_x(node_i.leaf_count);
  }
  
  hpart->root = n_tip + 1;
  
  return Rcpp::XPtr<TreeDist::HPart>(hpart, true);
}

#include <Rcpp.h>
#include "hpart.h"

using namespace Rcpp;

// Forward declaration
size_t build_node_from_list(const RObject& node,
                            std::vector<TreeDist::HNode>& nodes,
                            const int n_tip,
                            int& next_index,
                            const size_t n_block);

// [[Rcpp::export]]
SEXP build_hpart_from_list(RObject tree, const int n_tip) {
  const size_t vec_size = n_tip + n_tip + 2; // max nodes for a binary tree + 1
  
  auto hpart = new TreeDist::HPart();
  hpart->nodes.resize(vec_size);
  
  const size_t n_block = (n_tip + 63) / 64;
  const double n_tip_recip = 1 / static_cast<double>(n_tip);
  
  // Initialize leaves and internal nodes
  for (size_t i = 1; i < vec_size; ++i) {
    TreeDist::HNode& n = hpart->nodes[i];
    n.n_tip = n_tip;
    n.n_tip_reciprocal = n_tip_recip;
    n.bitset.assign(n_block, 0ULL);
    n.leaf_count = 0;
    n.all_kids_leaves = true;
    n.label = -1;
    n.children.clear();
  }
  
  // Start internal nodes at n_tip + 1
  int next_index = n_tip + 1;
  
  // Build recursively; returns index into nodes vector (1-based)
  hpart->root = build_node_from_list(tree, hpart->nodes, n_tip, next_index, n_block);
  
  return Rcpp::XPtr<TreeDist::HPart>(hpart, true);
}


// Recursive builder
size_t build_node_from_list(const RObject& node,
                            std::vector<TreeDist::HNode>& nodes,
                            int n_tip,
                            int& next_index,
                            size_t n_block) {
  
  if (Rf_isInteger(node) || Rf_isReal(node)) {
    const IntegerVector leaf_vec(node);
    if (leaf_vec.size() != 1) {
      Rcpp::stop("List must only contain integers, not vectors of integers");
    }
    const int leaf_label = leaf_vec[0];         // 1-based R leaf label
    const size_t leaf_idx = leaf_label - 1;     // 0-based label for HNode
    TreeDist::HNode& leaf = nodes[leaf_label];
    leaf.label = leaf_idx;
    leaf.leaf_count = 1;
    leaf.bitset[leaf_idx / 64] = 1ULL << (leaf_idx % 64);
    leaf.all_kids_leaves = true;
    return leaf_label;
  } else if (Rf_isVectorList(node)) {
    int n_elements = Rf_length(node);
    
    // Special case: a single leaf inside a length-1 list
    if (n_elements == 1 &&
        (Rf_isInteger(VECTOR_ELT(node,0)) || Rf_isReal(VECTOR_ELT(node,0)))) {
      return build_node_from_list(VECTOR_ELT(node, 0), nodes, n_tip, next_index, n_block);
    }
    
    // Allocate a new internal node
    const size_t my_idx = static_cast<size_t>(next_index++);
    TreeDist::HNode& n = nodes[my_idx];
    n.children.reserve(n_elements);
    n.leaf_count = 0;
    n.all_kids_leaves = true;
    n.n_tip = n_tip;
    
    // Recurse over children
    for (int i = 0; i < n_elements; ++i) {
      SEXP child = VECTOR_ELT(node, i);
      const size_t child_idx = build_node_from_list(child, nodes, n_tip, next_index, n_block);
      n.children.push_back(child_idx);
      
      // Merge bitsets
      std::transform(
        n.bitset.begin(), n.bitset.end(),
        nodes[child_idx].bitset.begin(),
        n.bitset.begin(),
        std::bit_or<uint64_t>()
      );
      
      n.leaf_count += nodes[child_idx].leaf_count;
      if (nodes[child_idx].leaf_count > 1) n.all_kids_leaves = false;
    }
    
    n.x_log_x = TreeDist::x_log_x(n.leaf_count);
    
    return my_idx;
  }
  
  // Invalid node type
  Rcpp::stop("Invalid node type");
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
    for (const size_t cidx : hp->nodes[i].children) {
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
