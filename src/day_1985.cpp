#include <Rcpp/Lightest>
using namespace Rcpp;

#include "tree_distances.h" /* includes <TreeTools/SplitList.h> */
#include "information.h"

#include <TreeTools.h> /* for root_on_node() */
#include <TreeTools/root_tree.h> /* for root_on_node() */
#include <TreeTools/ClusterTable.h> /* for ClusterTable() */
using TreeTools::ClusterTable;
using TreeTools::ct_stack_threshold;

#include <cmath> /* for log2(), ceil() */

struct StackEntry { int32 L, R, N, W; };

// COMCLUSTER computes a strict consensus tree in O(knn).
// COMCLUST requires O(kn).
// trees is a list of objects of class phylo.
// [[Rcpp::export]]
int COMCLUST(const List& trees) {
  
  int32 v = 0;
  int32 w = 0;
  int32 L, R, N, W;

  ClusterTable X(List(trees(0)));
  const int32 n_tip = X.N();

  StackEntry* S_ptr;
  std::array<StackEntry, ct_stack_threshold> S_stack;
  std::vector<StackEntry> S_heap;

  if (n_tip <= ct_stack_threshold) {
    S_ptr = S_stack.data();
  } else {
    S_heap.resize(n_tip);
    S_ptr = S_heap.data();
  }

  for (int32 i = 1; i < trees.length(); ++i) {
    int32 Spos = 0; // Empty the stack S
    
    X.CLEAR();
    ClusterTable Ti(List(trees(i)));
    Ti.TRESET();
    Ti.NVERTEX(&v, &w);
    
    do {
      if (Ti.is_leaf(v)) {
        S_ptr[Spos++] = {X.ENCODE(v), X.ENCODE(v), 1, 1};
      } else {
        const StackEntry& top = S_ptr[--Spos];
        L = top.L;
        R = top.R;
        N = top.N;
        const int32 W_i = top.W;
        
        W = 1 + W_i;
        w = w - W_i;
        while (w) {
          const StackEntry &entry = S_ptr[--Spos];
          L = std::min(L, entry.L); // Faster than ternary operator
          R = std::max(R, entry.R);
          N += entry.N;
          W += entry.W;
          w -= entry.W;
        };

        S_ptr[Spos++] = {L, R, N, W};

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

// trees is a list of objects of class phylo, all with the same tip labels
// (try RenumberTips(trees, trees[[1]]))
template <typename StackContainer>
double calc_consensus_info(const List &trees, const LogicalVector &phylo,
                           const NumericVector& p, StackContainer& S) {
  int32 v = 0;
  int32 w = 0;
  int32 L, R, N, W;

  const int32 n_trees = trees.length();

  std::vector<ClusterTable> tables;
  if (std::size_t(n_trees) > tables.max_size()) {
    Rcpp::stop("Not enough memory available to compute consensus of so many trees"); // LCOV_EXCL_LINE
  }

  tables.reserve(n_trees);
  for (int32 i = n_trees; i--; ) {
    tables.emplace_back(ClusterTable(List(trees(i))));
  }

  const int32 n_tip = tables[0].N();
  const int32 thresh = p[0] <= 0.5 ?
      (n_trees / 2) + 1 : // Splits must occur in MORE THAN 0.5 to be in majority.
      std::ceil(p[0] * n_trees);
  const int32 must_occur_before = 1 + n_trees - thresh;
  
  std::array<int32, ct_stack_threshold> split_count_stack;
  std::vector<int32> split_count_heap;
  int32* split_count;
  if (n_tip < ct_stack_threshold) {
    split_count = split_count_stack.data();
  } else {
    split_count_heap.resize(n_tip);
    split_count = split_count_heap.data();
  }

  StackEntry *const S_start = S.data();
  const bool phylo_info = phylo[0];
  double info = 0;
  
  const std::size_t ntip_3 = n_tip - 3;
  // All clades in p consensus must occur in first (1-p) of trees.
  for (int32 i = 0; i < must_occur_before; ++i) {
    if (tables[i].NOSWX(ntip_3)) {
      continue;
    }
    
    std::vector<int32> split_size(n_tip);
    std::fill(split_count, split_count + n_tip, 1);
    
    for (int32 j = i + 1; j < n_trees; ++j) {
      
      tables[i].CLEAR();
      
      tables[j].TRESET();
      tables[j].READT(&v, &w);
      
      int32 j_pos = 0;
      StackEntry* S_top = S_start; // Empty the stack S
      
      do {
        if (IS_LEAF(v)) {
          const auto enc_v = tables[i].ENCODE(v);
          *S_top++ = {enc_v, enc_v, 1, 1};
        } else {
          const StackEntry& entry = *--S_top;
          L = entry.L;
          R = entry.R;
          N = entry.N;
          W = 1 + entry.W;
          w -= entry.W;
          while (w) {
            const StackEntry& next = *--S_top;         
            L = std::min(L, next.L); // Faster than ternary operator
            R = std::max(R, next.R);
            N += next.N;
            W += next.W;
            w -= next.W;
          };
          *S_top++ = {L, R, N, W};
          
          ++j_pos;
          if (tables[j].GETSWX(&j_pos)) {
            // Split has already been counted; next!
          } else {
            if (N == R - L + 1) { // L..R is contiguous, and must be tested
              if (tables[i].CLUSTONL(L, R)) {
                tables[j].SETSWX(j_pos);
                assert(L > 0);
                ++split_count[L - 1];
                if (!split_size[L - 1]) {
                  split_size[L - 1] = N;
                }
                assert(split_size[L - 1] > 0);
              } else if (tables[i].CLUSTONR(L, R)) {
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
    
    int32 splits_found = 0;
    for (int32 k = n_tip; k--; ) {
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
IntegerVector robinson_foulds_all_pairs(const List& tables) {
  const int n_trees = static_cast<int>(tables.size());
  if (n_trees < 2) return IntegerVector(0);

  std::vector<ClusterTable*> tbl;
  tbl.reserve(n_trees);
  for (int i = 0; i < n_trees; ++i) {
    Rcpp::XPtr<ClusterTable> xp = tables[i];
    tbl.push_back(xp.get()); // .get() on XPtr => ClusterTable*
  }

  const size_t n_pairs = static_cast<size_t>(n_trees) * (n_trees - 1) / 2;
  IntegerVector shared = Rcpp::no_init(n_pairs);
  int *write_pos = INTEGER(shared); // direct pointer into R memory

  const int32 n_tip = tbl[0]->N();
  StackEntry* S_start;
  std::array<StackEntry, ct_stack_threshold> S_stack;
  std::vector<StackEntry> S_heap;
  if (n_tip <= ct_stack_threshold) {
    S_start = S_stack.data();
  } else {
    S_heap.resize(n_tip);
    S_start = S_heap.data();
  }

  for (int i = 0; i < n_trees - 1; ++i) {
    
    ClusterTable* Xi = tbl[i];
    
    for (int j = i + 1; j < n_trees; ++j) {
      
      int32 v;
      int32 w;
      int32 n_shared = 0;
      
      ClusterTable* Tj = tbl[j];
      
      // Reset stack pointer for each tree pair comparison
      StackEntry* S_top = S_start;
      
      Tj->TRESET();
      Tj->NVERTEX_short(&v, &w);
      
      while (v) {
        if (Tj->is_leaf(v)) {
          const auto enc_v = Xi->ENCODE(v);
          *S_top++ = {enc_v, enc_v, 1, 1};
        } else {
          const StackEntry& entry = *--S_top;
          int32 L = entry.L;
          int32 R = entry.R;
          int32 N = entry.N;
          const int32 W_i = entry.W;
          int32 W = 1 + W_i;
          
          w -= W_i;
          
          if (w) { // Unroll first iteration - common case
            const StackEntry& entry = *--S_top;
            const int32 W_i = entry.W;
            
            L = std::min(L, entry.L); // Faster than ternary operator
            R = std::max(R, entry.R);
            N += entry.N;
            W += W_i;
            w -= W_i;
            
            while (w) {
              const StackEntry& entry = *--S_top;
              const int32 W_i = entry.W;
              
              L = std::min(L, entry.L);
              R = std::max(R, entry.R);
              N += entry.N;
              W += W_i;
              w -= W_i;
            }
          }
          
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

// [[Rcpp::export]]
double consensus_info(const List trees, const LogicalVector phylo,
                      const NumericVector p) {
  if (p[0] > 1 + 1e-15) { // epsilon catches floating point error
    Rcpp::stop("p must be <= 1.0 in consensus_info()");
  } else if (p[0] < 0.5) {
    Rcpp::stop("p must be >= 0.5 in consensus_info()");
  }

  // First, peek at the tree size to determine allocation strategy
  // We'll create a temporary ClusterTable just to check the size
  try {
    ClusterTable temp_table(Rcpp::List(trees(0)));
    const int32 n_tip = temp_table.N();
    
    if (n_tip <= ct_stack_threshold) {
      // Small tree: use stack-allocated array
      std::array<StackEntry, ct_stack_threshold> S;
      return calc_consensus_info(trees, phylo, p, S);
    } else {
      // Large tree: use heap-allocated vector
      std::vector<StackEntry> S(n_tip);
      return calc_consensus_info(trees, phylo, p, S);
    }
  } catch(const std::exception& e) {
    Rcpp::stop(e.what());
  }
  
  ASSERT(false && "Unreachable code in consensus_tree");
  return 0.0;
}
