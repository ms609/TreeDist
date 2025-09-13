#include "hpart.h"
#include <Rcpp.h>
using namespace Rcpp;

namespace TreeDist {

// --- postorder recompute internal nodes ---
void recompute_bitsets_postorder(TreeDist::HPart* hpart,
                                 const size_t node_idx,
                                 const std::vector<int> &mapping) {
  auto &node = hpart->nodes[node_idx];
  uint64_t* node_bits = node.bitset;
  const size_t n_block = hpart->n_block;
  
  if (node.children.empty()) {
    // Leaf node
    if (node.leaf_count != 1) {
      Rcpp::stop("Leaf node has leaf_count != 1");
    }
    int new_index = mapping[node.label]; // mapping is 0-based
    node_bits[node.label / 64] = 0;
    
    node.label = new_index;
    const size_t block_idx = new_index >> 6; // = new_index / 64
    const size_t bit_idx = new_index & 63;   // = new_index % 64
    node_bits[block_idx] = 1ULL << bit_idx;
    
  } else {
    // Postorder: compute first child
    recompute_bitsets_postorder(hpart, node.children[0], mapping);
    auto &first_child = hpart->nodes[node.children[0]];
    const uint64_t* child_bits = first_child.bitset;
    
    // Copy first child’s bitset into parent
    std::copy(child_bits, child_bits + n_block, node_bits);
    
    node.leaf_count = first_child.leaf_count;
    node.all_kids_leaves = (node.leaf_count == 1);
    
    // Remaining children
    for (size_t ci = 1; ci < node.children.size(); ++ci) {
      size_t child_idx = node.children[ci];
      recompute_bitsets_postorder(hpart, child_idx, mapping);
      auto &child = hpart->nodes[child_idx];
      const uint64_t* cbits = child.bitset;
      
      // OR child into parent in-place
      std::transform(node_bits, node_bits + n_block,
                     cbits, node_bits,
                     std::bit_or<>{});
      
      node.leaf_count += child.leaf_count;
      if (child.leaf_count > 1) {
        node.all_kids_leaves = false;
      }
    }
  }
}

void relabel_hpart_xptr(TreeDist::HPart* hpart, const std::vector<int>& map) {
  recompute_bitsets_postorder(hpart, hpart->root, map);
}

} // namespace TreeDist


// [[Rcpp::export]]
void relabel_hpart(SEXP hpart_ptr, const std::vector<int>& map) {
  Rcpp::XPtr<TreeDist::HPart> hpart_xptr(hpart_ptr);
  TreeDist::HPart* hpart = hpart_xptr.get();
  TreeDist::relabel_hpart_xptr(hpart, map);
}

