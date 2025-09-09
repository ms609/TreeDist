#include "hpart.h"
#include <Rcpp.h>

// [[Rcpp::depends(TreeTools)]]
#include <TreeTools/renumber_tree.h> // for preorder_edges_and_nodes

using namespace Rcpp;
using Node = TreeDist::HNode;

// Forward
static void compute_bitsets_from(int node, HPart &hp);

// Recursive build
static void build_node(int node, HPart &hp,
                       const std::vector<std::vector<int>> &children,
                       int nTip) {
  Node &nd = hp.nodes[node];
  if (node <= nTip) {
    nd.labelIndex = node - 1; // 0-based
  } else {
    for (int c : children[node]) {
      build_node(c, hp, children, nTip);
      nd.children.push_back(c);
      if (c > nTip) {
        nd.allKidsLeaves = false;
      }
    }
  }
}

static void compute_bitsets_from(int node, HPart &hp) {
  Node &nd = hp.nodes[node];
  nd.bitset.assign(hp.nBlocks, 0);
  
  if (nd.labelIndex >= 0) {
    int idx = nd.labelIndex;
    nd.bitset[idx / 64] |= (1ULL << (idx % 64));
    nd.leafCount = 1;
  } else {
    for (int c : nd.children) {
      compute_bitsets_from(c, hp);
      Node &ch = hp.nodes[c];
      for (int b = 0; b < hp.nBlocks; ++b) {
        nd.bitset[b] |= ch.bitset[b];
      }
      nd.leafCount += ch.leafCount;
    }
  }
}

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
  }
  
  // Traverse nodes in postorder
  for (size_t i = vec_size; i > (size_t)n_tip; --i) {
    hpart->nodes[i].children.reserve(children[i].size());
    for (size_t child_id : children[i]) {
      hpart->nodes[i].children.push_back(&hpart->nodes[child_id]);
      for (size_t chunk = 0; chunk < nodes[i].bitset.size(); ++chunk) {
        hpart->nodes[i].bitset[chunk] |= hpart->nodes[child_id].bitset[chunk];
      }
    }
  }
  
  hpart->root = hpart->nodes[n_tip + 1];
  
  return Rcpp::XPtr<HPart>(hpart, true);
}
