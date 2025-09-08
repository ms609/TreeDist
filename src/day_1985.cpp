#include <Rcpp/Lightest>
using namespace Rcpp;

#include "tree_distances.h" /* includes <TreeTools/SplitList.h> */
#include "information.h"

#include <TreeTools.h> /* for root_on_node() */
#include <TreeTools/root_tree.h> /* for root_on_node() */
#include <TreeTools/ClusterTable.h> /* for ClusterTable() */
using TreeTools::ClusterTable;
using TreeTools::ct_max_leaves;

#include <array> /* for array */
#include <bitset> /* for bitset */
#include <vector> /* for vector */
#include <cmath> /* for log2(), ceil() */
#include <memory> /* for unique_ptr, make_unique */

struct StackEntry { int16 L, R, N, W; };

// COMCLUSTER computes a strict consensus tree in O(knn).
// COMCLUST requires O(kn).
// trees is a list of objects of class phylo.
// [[Rcpp::export]]
int COMCLUST(List trees) {
  
  int16 v = 0, w = 0;
  int16 L, R, N, W;
  int16 L_i, R_i, N_i, W_i;
  
  ClusterTable X(List(trees(0)));
  std::array<int16, TreeTools::ct_max_leaves> S;
  
  for (int16 i = 1; i != trees.length(); i++) {
    int16 Spos = 0; // Empty the stack S
    
    X.CLEAR();
    ClusterTable Ti(List(trees(i)));
    Ti.TRESET();
    Ti.NVERTEX(&v, &w);
    
    do {
      if (Ti.is_leaf(v)) {
        CT_PUSH(X.ENCODE(v), X.ENCODE(v), 1, 1);
      } else {
        CT_POP(L, R, N, W_i);
        W = 1 + W_i;
        w = w - W_i;
        while (w) {
          CT_POP(L_i, R_i, N_i, W_i);
          if (L_i < L) L = L_i;
          if (R_i > R) R = R_i;
          N += N_i;
          W += W_i;
          w -= W_i;
        };
        CT_PUSH(L, R, N, W);
        if (N == R - L + 1) { // L..R is contiguous, and must be tested
          X.SETSW(L, R);
        }
      }
      Ti.NVERTEX(&v, &w);
    } while (v);
    X.UPDATE();
  }
  
  return X.SHARED() - 2; // Subtract All-tips & All-ingroup
}

#define IS_LEAF(a) (a) <= n_tip
  
// COMCLUSTER computes a strict consensus tree in O(knn).
// COMCLUST requires O(kn).
// trees is a list of objects of class phylo, all with the same tip labels
// (try RenumberTips(trees, trees[[1]]))
// [[Rcpp::export]]
double consensus_info(const List trees, const LogicalVector phylo,
                      const NumericVector p) {
  
  int16 v = 0, w = 0,
    L, R, N, W,
    L_j, R_j, N_j, W_j
  ;
  const int16 n_trees = trees.length();
  
  std::vector<ClusterTable> tables;
  if (std::size_t(n_trees) > tables.max_size()) {
    Rcpp::stop("Not enough memory available to compute consensus of so many trees"); // LCOV_EXCL_LINE
  }
  tables.reserve(n_trees);
  for (int16 i = n_trees; i--; ) {
    tables.emplace_back(ClusterTable(List(trees(i))));
  }
  
  if (p[0] > 1) {
    Rcpp::stop("p must be <= 1.0 in consensus_info()");
  } else if (p[0] < 0.5) {
    Rcpp::stop("p must be >= 0.5 in consensus_info()");
  }
  const int16
    n_tip = tables[0].N(),
    thresh = p[0] <= 0.5 ?
      (n_trees / 2) + 1 : // Splits must occur in MORE THAN 0.5 to be in majority.
      std::ceil(p[0] * n_trees),
    must_occur_before = 1 + n_trees - thresh
  ;
  
  const bool phylo_info = phylo[0];
  
  std::array<int16, TreeTools::ct_stack_size * TreeTools::ct_max_leaves> S;
  std::array<int16, TreeTools::ct_max_leaves> split_count;
  
  double info = 0;
  
  const std::size_t ntip_3 = n_tip - 3;
  // All clades in p consensus must occur in first (1-p) of trees.
  for (int16 i = 0; i != must_occur_before; i++) {
    if (tables[i].NOSWX(ntip_3)) {
      continue; 
    }
    
    std::vector<int16> split_size(n_tip);
    std::fill(split_count.begin(), split_count.begin() + n_tip, 1);
    
    for (int16 j = i + 1; j != n_trees; j++) {
      
      tables[i].CLEAR();
      
      tables[j].TRESET();
      tables[j].READT(&v, &w);
      
      int16 j_pos = 0, Spos = 0; // Empty the stack S
      
      do {
        if (IS_LEAF(v)) {
          CT_PUSH(tables[i].ENCODE(v), tables[i].ENCODE(v), 1, 1);
        } else {
          CT_POP(L, R, N, W_j);
          W = 1 + W_j;
          w = w - W_j;
          while (w) {
            CT_POP(L_j, R_j, N_j, W_j);
            if (L_j < L) L = L_j;
            if (R_j > R) R = R_j;
            N = N + N_j;
            W = W + W_j;
            w = w - W_j;
          };
          CT_PUSH(L, R, N, W);
          
          ++j_pos;
          if (tables[j].GETSWX(&j_pos)) {
            // Split has already been counted; next!
          } else {
            if (N == R - L + 1) { // L..R is contiguous, and must be tested
              if (tables[i].CLUSTONL(&L, &R)) {
                tables[j].SETSWX(j_pos);
                assert(L > 0);
                ++split_count[L - 1];
                if (!split_size[L - 1]) {
                  split_size[L - 1] = N;
                }
                assert(split_size[L - 1] > 0);
              } else if (tables[i].CLUSTONR(&L, &R)) {
                tables[j].SETSWX(j_pos);
                assert(R > 0);
                ++split_count[R - 1];
                if (!split_size[R - 1]) {
                  split_size[R - 1] = N;
                }
                assert(split_size[R - 1] > 0);
              }
            }
          }
        }
        tables[j].NVERTEX_short(&v, &w);
      } while (v);
    }
    
    int16 splits_found = 0;
    for (int16 k = n_tip; k--; ) {
      if (split_count[k] >= thresh) {
        ++splits_found;
        if (phylo_info) {
          info += split_phylo_info(split_size[k], &n_tip, 
                                   split_count[k] / double(n_trees));
        } else {
          info += split_clust_info(split_size[k], &n_tip, 
                                   split_count[k] / double(n_trees));
        }
        // If we have a perfectly resolved tree, break.
        if (splits_found == n_tip - 3) {
          return phylo_info ? info : info * n_tip;
        }
      }
    }
  }
  
  // Convert clustering entropy to *total* information 
  return phylo_info ? info : info * n_tip;
}

// [[Rcpp::export]]
IntegerVector robinson_foulds_all_pairs(List tables) {
  const int n_trees = static_cast<int>(tables.size());
  if (n_trees < 2) return IntegerVector(0);
  
  std::vector<Rcpp::XPtr<ClusterTable>> xptrs;
  xptrs.reserve(n_trees);
  for (int i = 0; i < n_trees; ++i) {
    Rcpp::XPtr<ClusterTable> xp = tables[i];
    xptrs.emplace_back(xp);
  }
  
  std::vector<ClusterTable*> tbl;
  tbl.reserve(n_trees);
  for (int i = 0; i < n_trees; ++i) {
    tbl.push_back(xptrs[i].get()); // .get() on XPtr => ClusterTable*
  }
  
  const size_t n_pairs = static_cast<size_t>(n_trees) * (n_trees - 1) / 2;
  IntegerVector shared = Rcpp::no_init(n_pairs);
  int *write_pos = INTEGER(shared); // direct pointer into R memory
  
  std::array<StackEntry, ct_max_leaves> S_entries;
  StackEntry* S_top = S_entries.data(); // stack top pointer (points one past last element)
  
  for (int i = 0; i < n_trees - 1; ++i) {
    
    ClusterTable* Xi = tbl[i];
    
    for (int j = i + 1; j < n_trees; ++j) {
      
      int16 v, w;
      int16 n_shared = 0;
      
      ClusterTable* Tj = tbl[j];
      
      // Reset stack pointer for each tree pair comparison
      S_top = S_entries.data();
      
      Tj->TRESET();
      Tj->NVERTEX_short(&v, &w);
      
      while (v) {
        if (Tj->is_leaf(v)) {
          const auto enc_v = Xi->ENCODE(v);
          ASSERT(S_top < S_entries.data() + ct_max_leaves);
          *S_top++ = {enc_v, enc_v, 1, 1};
        } else {
          ASSERT(S_top > S_entries.data());
          const StackEntry& entry = *--S_top;
          int16 L = entry.L;
          int16 R = entry.R;
          int16 N = entry.N;
          const int16 W_i = entry.W;
          int16 W = 1 + W_i;
          
          w -= W_i;
          
          if (w) { // Unroll first iteration - common case
            ASSERT(S_top > S_entries.data());
            const StackEntry& entry = *--S_top;
            const int16 W_i = entry.W;
            
            L = std::min(L, entry.L); // Faster than ternary operator
            R = std::max(R, entry.R);
            N += entry.N;
            W += W_i;
            w -= W_i;
            
            while (w) {
              ASSERT(S_top > S_entries.data());
              const StackEntry& entry = *--S_top;
              const int16 W_i = entry.W;
              
              L = std::min(L, entry.L);
              R = std::max(R, entry.R);
              N += entry.N;
              W += W_i;
              w -= W_i;
            }
          }
          
          ASSERT(S_top < S_entries.data() + ct_max_leaves);
          *S_top++ = {L, R, N, W};
          
          if (N == R - L + 1) { // L..R is contiguous, and must be tested
            if (Xi->ISCLUST(L, R)) ++n_shared;
          }
        }
        Tj->NVERTEX_short(&v, &w); // Doesn't count all-ingroup or all-tips
      }
      *write_pos++ = n_shared;
    }
  }
  
  return shared;
}
