#include <Rcpp.h>
using namespace Rcpp;

#include "tree_distances.h" /* includes <TreeTools/SplitList.h> */
#include "information.h"

#include <TreeTools.h> /* for root_on_node() */
#include <TreeTools/root_tree.h> /* for root_on_node() */
#include <TreeTools/ClusterTable.h> /* for ClusterTable() */
using TreeTools::ClusterTable;

#include <array> /* for array */
#include <bitset> /* for bitset */
#include <vector> /* for vector */
#include <cmath> /* for log2(), ceil() */
#include <memory> /* for unique_ptr, make_unique */

// Modelled on https://CRAN.R-project.org/package=Rcpp/vignettes/Rcpp-modules.pdf
// [[Rcpp::export]]
SEXP ClusterTable_new(List phylo) {
  XPtr<ClusterTable> ptr(new ClusterTable (phylo), true);
  return ptr;
}

// [[Rcpp::export]]
IntegerMatrix ClusterTable_matrix(SEXP xp) {
  XPtr<ClusterTable> ptr(xp);
  return ptr->X_contents();
}

// [[Rcpp::export]]
IntegerVector ClusterTable_decode(SEXP xp) {
  XPtr<ClusterTable> ptr(xp);
  return ptr->X_decode();
}

// COMCLUSTER computes a strict consensus tree in O(knn).
// COMCLUST requires O(kn).
// trees is a list of objects of class phylo.
// [[Rcpp::export]]
int COMCLUST (List trees) {
  
  int16 v = 0, w = 0,
    L, R, N, W,
    L_i, R_i, N_i, W_i
  ;
  
  ClusterTable X(List(trees(0)));
  std::array<int16, CT_MAX_LEAVES> S;
  
  for (int16 i = 1; i != trees.length(); i++) {
    int16 Spos = 0; // Empty the stack S
    
    X.CLEAR();
    ClusterTable Ti(List(trees(i)));
    Ti.TRESET();
    Ti.NVERTEX(&v, &w);
    
    do {
      if (Ti.is_leaf(&v)) {
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
          X.SETSW(&L, &R);
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
double consensus_info (const List trees, const LogicalVector phylo,
                       const NumericVector p) {
  
  int16 v = 0, w = 0,
    L, R, N, W,
    L_j, R_j, N_j, W_j
  ;
  const int16 n_trees = trees.length();
  
  std::vector<ClusterTable> tables;
  tables.reserve(n_trees);
  for (int16 i = n_trees; i--; ) {
    tables.emplace_back(ClusterTable(List(trees(i))));
  }
  
  if (p[0] > 1) {
    throw std::range_error("p must be <= 1.0 in consensus_info()");
  } else if (p[0] < 0.5) {
    throw std::range_error("p must be >= 0.5 in consensus_info()");
  }
  const int16
    n_tip = tables[0].N(),
    thresh = p[0] <= 0.5 ?
      (n_trees / 2) + 1 : // Splits must occur in MORE THAN 0.5 to be in majority.
      std::ceil(p[0] * n_trees),
    must_occur_before = 1 + n_trees - thresh
  ;
  
  const bool phylo_info = phylo[0];
  
  std::array<int16, CT_STACK_SIZE * CT_MAX_LEAVES> S;
  std::array<int16, CT_MAX_LEAVES> split_count;
  
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
                tables[j].SETSWX(&j_pos);
                assert(L > 0);
                ++split_count[L - 1];
                if (!split_size[L - 1]) {
                  split_size[L - 1] = N;
                }
                assert(split_size[L - 1] > 0);
              } else if (tables[i].CLUSTONR(&L, &R)) {
                tables[j].SETSWX(&j_pos);
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
  int16 
    v = 0, w = 0,
    L, R, N, W,
    L_i, R_i, N_i, W_i
  ;
  const int16 n_trees = tables.length();
  if (n_trees < 2) return IntegerVector(0);
  
  IntegerVector shared(n_trees * (n_trees - 1) / 2);
  IntegerVector::iterator write_pos = shared.begin();
  std::array<int16, CT_MAX_LEAVES> S;
  
  for (int16 i = 0; i != n_trees - 1; i++) {
    Rcpp::XPtr<ClusterTable> table_i = tables(i);
    Rcpp::XPtr<ClusterTable> Xi(table_i);
    
    for (int16 j = i + 1; j != n_trees; j++) {
      Rcpp::XPtr<ClusterTable> table_j = tables(j);
      Rcpp::XPtr<ClusterTable> Tj(table_j);
      int16 Spos = 0; // Empty the stack S
      int16 n_shared = 0;
      
      Tj->TRESET();
      Tj->NVERTEX_short(&v, &w);
      
      do {
        if (Tj->is_leaf(&v)) {
          CT_PUSH(Xi->ENCODE(v), Xi->ENCODE(v), 1, 1);
        } else {
          CT_POP(L, R, N, W_i);
          W = 1 + W_i;
          w = w - W_i;
          while (w) {
            CT_POP(L_i, R_i, N_i, W_i);
            if (L_i < L) L = L_i;
            if (R_i > R) R = R_i;
            N = N + N_i;
            W = W + W_i;
            w = w - W_i;
          };
          CT_PUSH(L, R, N, W);
          if (N == R - L + 1) { // L..R is contiguous, and must be tested
            if (Xi->ISCLUST(&L, &R)) ++n_shared;
          }
        }
        Tj->NVERTEX_short(&v, &w); // Doesn't count all-ingroup or all-tips
      } while (v);
      *write_pos++ = n_shared;
    }
  }
  return shared;
}
