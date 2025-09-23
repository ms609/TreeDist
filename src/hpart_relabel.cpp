#include "hpart.h"
#include <Rcpp.h>
using namespace Rcpp;

namespace TreeDist {

// --- postorder recompute internal nodes ---
void recompute_bitsets_postorder(TreeDist::HPart &hpart, const size_t node_idx,
                                 const std::vector<int> &mapping,
                                 const size_t n_block) {
  auto &node = hpart.nodes[node_idx];
  
  if (node.children.empty()) {
    // Leaf node
    if (node.leaf_count != 1) {
      Rcpp::stop("Leaf node has leaf_count != 1");                              // #nocov
    }
    int new_index = mapping[node.label]; // mapping is 0-based
    node.label = new_index;
    node.bitset.assign(n_block, 0ULL);
    node.bitset[new_index / 64] = 1ULL << (new_index % 64);
  } else {
    // Postorder: first child
    recompute_bitsets_postorder(hpart, node.children[0], mapping, n_block);
    auto &first_child = hpart.nodes[node.children[0]];
    node.bitset = first_child.bitset;
    node.leaf_count = first_child.leaf_count;
    node.all_kids_leaves = (node.leaf_count == 1);
    
    // Remaining children
    for (size_t ci = 1; ci < node.children.size(); ++ci) {
      size_t child_idx = node.children[ci];
      recompute_bitsets_postorder(hpart, child_idx, mapping, n_block);
      auto &child = hpart.nodes[child_idx];
      
      for (size_t chunk = 0; chunk < node.bitset.size(); ++chunk) {
        node.bitset[chunk] |= child.bitset[chunk];
      }
      
      node.leaf_count += child.leaf_count;
      if (child.leaf_count > 1) {
        node.all_kids_leaves = false;
      }
    }
  }
}

} // namespace TreeDist


// [[Rcpp::export]]
void relabel_hpart(SEXP hpart_ptr, const std::vector<int>& map) {
  Rcpp::XPtr<TreeDist::HPart> hpart_xptr(hpart_ptr);
  TreeDist::HPart &hpart = *hpart_xptr;
  
  const size_t n_block = hpart.nodes[1].bitset.size();
  
  TreeDist::recompute_bitsets_postorder(hpart, hpart.root, map, n_block);
}
