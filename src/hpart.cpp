#include "hpart.h"
#include <Rcpp.h>

// [[Rcpp::depends(TreeTools)]]
#include <TreeTools/renumber_tree.h> // for preorder_edges_and_nodes

using namespace Rcpp;

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
  int nTip = tip_label.size();
  int Nnode = phy["Nnode"];
  
  // preorder edges
  IntegerMatrix reordered = TreeTools::preorder_edges_and_nodes(edge(_,0), edge(_,1));
  
  std::vector<std::vector<int>> children(nTip + Nnode + 1);
  for (int i = 0; i < reordered.nrow(); ++i) {
    int p = reordered(i, 0);
    int c = reordered(i, 1);
    children[p].push_back(c);
  }
  
  HPart* hp = new HPart();
  hp->nTips = nTip;
  hp->nBlocks = (nTip + 63) / 64;
  hp->rootIndex = nTip + 1;
  hp->nodes.resize(nTip + Nnode + 1);
  
  build_node(hp->rootIndex, *hp, children, nTip);
  compute_bitsets_from(hp->rootIndex, *hp);
  
  return Rcpp::XPtr<HPart>(hp, true); // managed by R
}
