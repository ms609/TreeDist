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
  std::vector<size_t> child_store;
  const size_t n_edge = n_tip + n_node - 1;
  child_store.reserve(n_edge);
  
  std::vector<size_t> first_child(vec_size + 1, 0);
  
  for (size_t i = 0; i < (size_t)reordered.nrow(); ++i) {
    first_child[reordered(i, 0)]++;
  }
  
  size_t total = 0;
  for (size_t i = 0; i < vec_size; ++i) {
    size_t i_children = first_child[i];
    first_child[i] = total;
    total += i_children;
  }
  first_child[vec_size] = total;
  
  child_store.resize(total);
    
  std::vector<size_t> offset = first_child;
  for (size_t i = 0; i < (size_t)reordered.nrow(); ++i) {
    size_t p = reordered(i, 0);
    size_t c = reordered(i, 1);
    child_store[offset[p]++] = c;
  }
  
  TreeDist::HPart* hpart = new TreeDist::HPart(n_tip, vec_size);
  const int n_block = hpart->n_block;
  
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
    
    node_i.children = &child_store[first_child[i]];
    size_t start = first_child[i];
    size_t end = first_child[i + 1];
    node_i.n_children = end - start;
    
    
    // Reset node properties
    node_i.all_kids_leaves = true;
    node_i.leaf_count = 0;
    
    // Process each child
    for (size_t j = start; j < end; ++j) {
      size_t child_id = child_store[j];
      const auto &child_node = hpart->nodes[child_id];
      
      // Accumulate leaf count
      node_i.leaf_count += child_node.leaf_count;
      
      // Check if any child is internal
      if (child_node.leaf_count > 1) {
        node_i.all_kids_leaves = false;
      }
      
      // Combine bitsets
      std::transform(
        node_i.bitset, node_i.bitset + n_block,
        child_node.bitset,
        node_i.bitset,
        std::bit_or<uint64_t>()
      );
    }
    
    // Precompute x*log(x)
    node_i.x_log_x = TreeDist::x_log_x(node_i.leaf_count);
  }
  hpart->root = n_tip + 1;
  
  return Rcpp::XPtr<TreeDist::HPart>(hpart, true);
}

// Forward declaration
size_t build_node_from_list(const RObject& node,
                            std::vector<TreeDist::HNode>& nodes,
                            const int n_tip,
                            int& next_index,
                            std::vector<size_t>& child_store);

// [[Rcpp::export]]
SEXP build_hpart_from_list(RObject tree, const int n_tip) {
  // Maximum nodes: for a binary tree, max internal nodes = n_tip - 1, 
  // so vec_size = n_tip + n_tip + 1.  We add one because there is no node
  // numbered zero.
  const size_t vec_size = n_tip + n_tip + 2;
  
  TreeDist::HPart* hpart = new TreeDist::HPart(n_tip, vec_size);
  
  // Initialize leaves and internal nodes
  for (size_t i = 1; i < vec_size; ++i) {
    TreeDist::HNode& n = hpart->nodes[i];
    n.leaf_count = 0;
    n.all_kids_leaves = true;
    n.label = -1;
    n.n_children = 0;
    n.children = nullptr;
  }
  
  std::vector<size_t> child_store;    // will hold all children contiguously
  child_store.reserve(vec_size - 1);  // rough upper bound
  
  int next_index = n_tip + 1;  // internal nodes start here
  
  // Recursive build
  hpart->root = build_node_from_list(tree, hpart->nodes, n_tip, next_index,
                                     child_store);
  
  // Move child_store into HPart
  hpart->child_store = std::move(child_store);
  
  return Rcpp::XPtr<TreeDist::HPart>(hpart, true);
}

// Recursive builder
size_t build_node_from_list(const RObject& node,
                            std::vector<TreeDist::HNode>& nodes,
                            const int n_tip,
                            int& next_index,
                            std::vector<size_t>& child_store) {
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
    leaf.n_children = 0;
    leaf.children = nullptr;
    return leaf_label;
    
  } else if (Rf_isVectorList(node)) {
    int n_elements = Rf_length(node);
    
    // Special case: single leaf in a length-1 list
    if (n_elements == 1 &&
        (Rf_isInteger(VECTOR_ELT(node,0)) || Rf_isReal(VECTOR_ELT(node,0)))) {
      return build_node_from_list(VECTOR_ELT(node, 0), nodes, n_tip, next_index,
                                  child_store);
    }
    
    // Allocate a new internal node
    const size_t my_idx = next_index++;
    TreeDist::HNode& n = nodes[my_idx];
    n.leaf_count = 0;
    n.all_kids_leaves = true;
    n.n_tip = n_tip;
    
    size_t start_pos = child_store.size(); // Doesn't work; we need a first
                                  // pass to populate properly - follow the phylo
                                  // approach below.
    
    for (int i = 0; i < n_elements; ++i) {
      SEXP child = VECTOR_ELT(node, i);
      size_t child_idx = build_node_from_list(child, nodes, n_tip, next_index, child_store);
      child_store.push_back(child_idx);
      
      // Merge bitsets
      std::transform(
        n.bitset, n.bitset + n.n_block,
        nodes[child_idx].bitset,
        n.bitset,
        std::bit_or<uint64_t>()
      );
      
      n.leaf_count += nodes[child_idx].leaf_count;
      if (nodes[child_idx].leaf_count > 1) n.all_kids_leaves = false;
    }
    
    Rcpp::Rcout << "start pos: " << start_pos <<
    n.children = &child_store[start_pos];
    n.n_children = child_store.size() - start_pos;
    n.x_log_x = TreeDist::x_log_x(n.leaf_count);
    
    return my_idx;
  }
  
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
    const int parent_id = node_ids[i];
    const auto& node_i = hp->nodes[i];
    for (size_t j = 0; j < node_i.n_children; ++j) {
      const size_t cidx = node_i.children[j];
      const int child_id = node_ids[cidx];
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
