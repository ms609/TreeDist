#include <Rcpp.h>
using namespace Rcpp;
#include <TreeTools.h> /* for root_on_node() */
#include <TreeTools/root_tree.h> /* for root_on_node() */
#include <cmath> /* for log2() */
#include <memory> /* for unique_ptr, make_unique */
#include "tree_distances.h"
#include "SplitList.h"

const int16
  INF = INTX_MAX,
  UNINIT = -999
;

const int_fast32_t
  DAY_MAX_LEAVES = 16383,
  FACT_MAX = DAY_MAX_LEAVES + DAY_MAX_LEAVES + 5 + 1
;

const double log_2 = log(2);
double ldfact[FACT_MAX];
double l2rooted[DAY_MAX_LEAVES];
double *l2unrooted = &l2rooted[0] - 1;
double l2[DAY_MAX_LEAVES];
__attribute__((constructor))
  void compute_double_factorials() {
    ldfact[0] = 0;
    ldfact[1] = 0;
    l2rooted[0] = 0;
    l2rooted[1] = 0;
    l2rooted[2] = 0;
    l2[1] = 0;
    l2[2] = 1;
    
    for (int_fast32_t i = 2; i != FACT_MAX; i++) {
      ldfact[i] = ldfact[i - 2] + log2(i);
    }

    for (int_fast32_t i = 3; i != DAY_MAX_LEAVES; i++) {
      l2[i] = log2(i);
      l2rooted[i] = ldfact[i + i - 3];
      assert(l2unrooted[i] == ldfact[i + i - 5]);
    }
  }

inline double split_phylo_info (const int16 n_in, const int16 *n_tip,
                                const double p) {
  const int16 n_out = *n_tip - n_in;
  assert(p > 0);
  assert(p <= 1);
  assert(n_in > 1);
  assert(n_out > 1);
  if (p == 1) {
    return (l2unrooted[*n_tip] - l2rooted[n_in] - l2rooted[n_out]);
  } else {
    const double 
      q = 1 - p,
      l2n = l2unrooted[*n_tip],
      l2n_consistent = l2rooted[n_in] + l2rooted[n_out],
      l2p_consistent = l2n_consistent - l2n,
      l2p_inconsistent = log2(-expm1(l2p_consistent * log_2)),
      l2n_inconsistent = l2p_inconsistent + l2n
    ;
    
    return(l2n + 
           p * (log2(p) - l2n_consistent) +
           q * (log2(q) - l2n_inconsistent));
  }
}
  
inline double split_clust_info (const int16 n_in, const int16 *n_tip,
                                const double p) {
  const int16 n_out = *n_tip - n_in;
  assert(p > 0);
  assert(p <= 1);
  assert(n_in > 1);
  assert(n_out > 1);
  return -p * (
      (((n_in * l2[n_in]) + (n_out * l2[n_out])) / *n_tip) - l2[*n_tip]
    );
}


class ClusterTable {
  
  const int16
    L_COL = 0,
    R_COL = 1,
    SWITCH_COL = 2,
    X_COLS = 3
  ;
  int16
    n_edge,
    n_internal,
    n_leaves,
    
    n_shared = 0,
    enumeration = 0,
    v_j,
    Tlen,
    Tpos = 0,
    X_ROWS
  ;
  std::unique_ptr<int16[]> leftmost_leaf, T, Xarr, visited_nth, internal_label;
  
  public:
    ClusterTable(List); // i.e. PREPARE(T)
    
    inline bool is_leaf(const int16 *v) {
      return *v <= n_leaves;
    }
    
    inline const int16 edges() {
      return n_edge;
    }
    
    inline const int16 leaves() {
      return n_leaves;
    }
    
    inline void ENTER(int16 v, int16 w) {
      T.get()[Tpos++] = v;
      T.get()[Tpos++] = w;
    }
    
    inline void READT(int16 *v, int16 *w) {
      *v = T[Tpos++];
      *w = T[Tpos++];
    }
    
    inline int16 N() {
      return n_leaves;
    }
    
    inline int16 M() {
      return n_internal;
    }
    
    inline void TRESET() {
      // This procedure prepares T for an enumeration of its entries, 
      // beginning with the first entry. 
      Tpos = 0;
    }
    
    inline void NVERTEX(int16 *v, int16 *w) {
      if (Tpos != Tlen) {
        READT(v, w);
        v_j = *v;
      } else {
        *v = 0;
        *w = 0;
      }
    }
    
    inline void NVERTEX_short(int16 *v, int16 *w) {
      // Don't count all-tips or all-ingroup: vertices 0, ROOT, Ingp.
      if (Tpos != Tlen - (2 * 3)) {
        READT(v, w);
        v_j = *v;
      } else {
        *v = 0;
        *w = 0;
      }
    }
    
    inline int16 LEFTLEAF() {
      // If NVERTEX has returned entry <vj, wj> in T, the leftmost leaf in the
      // subtree rooted at vj has entry <vk, wk> where k = j - wj.
      // This function procedure returns Vk as its value.
      return leftmost_leaf[v_j - 1];
    }
    
    inline void SET_LEFTMOST(int16 index, int16 val) {
      leftmost_leaf[index - 1] = val;
    }
    
    inline int16 GET_LEFTMOST(int16 index) {
      return leftmost_leaf[index - 1];
    }
    
    // Procedures to manipulate cluster tables, per Table 4 of Day 1985.
    
    inline int16 ENCODE(const int16 v) {
      // This function procedure returns as its value the internal label 
      // assigned to leaf v
      // MS note: input = v; output = X[v, 3]
      return internal_label[v - 1];
    }
    
    inline int16 DECODE(const int16 internal_relabeling) {
      // MS: input = X[v, 3], output = v
      return visited_nth[internal_relabeling - 1];
    }
    
    inline void VISIT_LEAF (const int16* leaf, int16* n_visited) {
      visited_nth[(*n_visited)++] = *leaf;
      internal_label[*leaf - 1] = *n_visited;
    }
    
    IntegerVector X_decode() {
      IntegerVector ret(N());
      for (int16 i = n_leaves; i--; ) {
        ret(i) = DECODE(i + 1);
      }
      return ret;
    }
    
    inline int16 X(int16 row, int16 col) {
      if (row < 1) throw std::range_error("Trying to read before start of X");
      if (row > X_ROWS) {
        Rcout << "Row " << row << " > number of rows, " << X_ROWS << ".\n";
        throw std::range_error("Trying to read past end of X");
      }
      return Xarr[(row - 1) * X_COLS + col];
    }
    
    inline void setX(int16 row, int16 col, int16 value) {
      if (row < 1) throw std::range_error("Trying to write before start of X");
      if (row > X_ROWS) throw std::range_error("Trying to write past end of X");
      Xarr[(row - 1) * X_COLS + col] = value;
    }
    
    IntegerMatrix X_contents() {
      IntegerMatrix ret(X_ROWS, 2);
      for (int16 i = X_ROWS; i--; ) {
        ret(i, 0) = X(i + 1, L_COL);
        ret(i, 1) = X(i + 1, R_COL);
      }
      return ret;
    }
    
    inline bool CLUSTONL(int16* L, int16* R) {
      return X(*L, L_COL) == *L && X(*L, R_COL) == *R;
    }
    
    inline bool CLUSTONR(int16* L, int16* R) {
      return X(*R, L_COL) == *L && X(*R, R_COL) == *R;
    }
    
    inline bool ISCLUST(int16* L, int16* R) {
      // This function procedure returns value true if cluster <L,R> is in X;
      // otherwise it returns value false
      return CLUSTONL(L, R) || CLUSTONR(L, R);
    }
    
    inline void CLEAR() {
      // Each cluster in X has an associated switch that is either cleared or 
      // set. 
      // This procedure clears every cluster switch in X. 
      for (int16 i = X_ROWS; i--; ) {
        Xarr[i * X_COLS + SWITCH_COL] = 0;
      }
    }
    
    inline void SETSWX(int16* row) {
      setX(*row, SWITCH_COL, 1);
    }
    
    inline bool GETSWX(int16* row) {
      return X(*row, SWITCH_COL);
    }
    
    inline void SETSW(int16* L, int16* R) {
      // If <L,R> is a cluster in X, this procedure sets the cluster switch for <L,R>. 
      if (CLUSTONL(L, R)) {
        ++n_shared;
        SETSWX(L);
      } else if (CLUSTONR(L, R)) {
        ++n_shared;
        SETSWX(R);
      }
    }
    
    inline void UPDATE(){
      // This procadure inspects every cluster switch in X.
      // If the switch for cluster <L,R> is cleared, UPDATE deletes <L,R> 
      // from X; thereafter ISCLUST(X,L,R) will return the value false. 
      for (int16 i = X_ROWS; i--; ) {
        int16 ptr = (i * X_COLS);
        if (!(Xarr[ptr + SWITCH_COL])) {
          Xarr[ptr + L_COL] = -Xarr[ptr + L_COL]; // 0
          Xarr[ptr + R_COL] = -Xarr[ptr + R_COL]; // 0
        }
      }
    }
    
    inline int16 SHARED() {
      return n_shared;
    }
    
    inline void ADDSHARED() {
      ++n_shared;
    }
    
    inline void XRESET() {
      // This procedure prepares X for an enumeration of its clusters
      enumeration = 0;
    }
    
    inline void NCLUS(int16* L, int16* R) {
      // This procedure returns the next cluster <L,R> in the current 
      // enumeration of clusters in X.
      // If m clusters are in X, they are returned by the first m invocations 
      // of NCLUS after initialization by XRESET; thereafter NCLUS returns the
      // invalid cluster <0,0>. 
      *L = X(enumeration, 0);
      *R = X(enumeration, 1);
      enumeration++;
    }
    
};

ClusterTable::ClusterTable(List phylo) { 
  
  const List rooted = TreeTools::root_on_node(phylo, 1);
  const IntegerMatrix edge = rooted["edge"];
  
  // BEGIN
  n_internal = rooted["Nnode"]; // = M
  CharacterVector leaf_labels = rooted["tip.label"];
  if (leaf_labels.length() > DAY_MAX_LEAVES) {
    throw std::length_error("Tree has too many leaves. "
                            "Contact the 'TreeDist' maintainer.");
  }
  n_leaves = leaf_labels.length(); // = N
  n_edge = edge.nrow();
  Tlen = M() + N() + M() + N();
  T = std::make_unique<int16[]>(Tlen);
  
  leftmost_leaf = std::make_unique<int16[]>(N() + M());
  visited_nth = std::make_unique<int16[]>(n_leaves);
  internal_label = std::make_unique<int16[]>(n_leaves);
  int16 n_visited = 0;
  std::unique_ptr<int16[]> weights = std::make_unique<int16[]>(N() + M() + 1);
  
  for (int16 i = 1; i != n_leaves + 1; i++) {
    SET_LEFTMOST(i, i);
    weights[i] = 0;
  }
  for (int16 i = n_leaves + 1; i != N() + M() + 1; i++) {
    SET_LEFTMOST(i, 0);
    weights[i] = 0;
  }
  for (int16 i = n_edge; i--; ) {
    const int16
      parent_i = edge(i, 0),
      child_i = edge(i, 1)
    ;
    if (!GET_LEFTMOST(parent_i)) SET_LEFTMOST(parent_i, GET_LEFTMOST(child_i));
    if (is_leaf(&child_i)) {
      VISIT_LEAF(&child_i, &n_visited);
      weights[parent_i]++;
      ENTER(child_i, 0);
    } else {
      weights[parent_i] += 1 + weights[child_i];
      ENTER(child_i, weights[child_i]);
    }
  }
  ENTER(edge(0, 0), weights[edge(0, 0)]);
  
  // BUILD Cluster table
  X_ROWS = n_leaves;
  Xarr = std::make_unique<int16[]>(X_COLS * X_ROWS);
  
  // This procedure constructs in X descriptions of the clusters in a
  // rooted tree described by the postorder sequence T with weights,
  // BUILD assigns each leaf an internal label so that every cluster
  // is a set {i : L ~ i ~ R] of internal labels; thus each cluster is
  // simply described by a pair <L,R> of internal labels.
  
  TRESET();
  for (int16 i = 1; i != N(); i++) {
    setX(i, L_COL, 0);
    setX(i, R_COL, 0);
  }
  int16 leafcode = 0, v, w, L, R = UNINIT, loc;
  
  NVERTEX(&v, &w);
  while (v) {
    if (is_leaf(&v)) {
      leafcode++;
      // We prepared the encoder in an earlier step, so need no X[v, 3] <- leafcode
      R = leafcode;
      NVERTEX(&v, &w);
    } else {
      L = ENCODE(LEFTLEAF());
      NVERTEX(&v, &w);
      loc = w == 0 ? R : L;
      setX(loc, L_COL, L);
      setX(loc, R_COL, R);
    }
  }
}

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

inline void push (int16 a, int16 b, int16 c, int16 d, std::unique_ptr<int16[]> &S, int16* Spos) {
  S.get()[(*Spos)++] = a;
  S.get()[(*Spos)++] = b;
  S.get()[(*Spos)++] = c;
  S.get()[(*Spos)++] = d;
}

inline void push (int16 a, int16 b, int16 c, int16 d, IntegerVector &S, int16* Spos) {
  S[(*Spos)++] = a;
  S[(*Spos)++] = b;
  S[(*Spos)++] = c;
  S[(*Spos)++] = d;
}

inline void pop (int16 *a, int16 *b, int16 *c, int16 *d, std::unique_ptr<int16[]> &S, int16* Spos) {
  *d = S.get()[--(*Spos)];
  *c = S.get()[--(*Spos)];
  *b = S.get()[--(*Spos)];
  *a = S.get()[--(*Spos)];
}

inline void pop (int16 *a, int16 *b, int16 *c, int16 *d, IntegerVector &S, int16* Spos) {
  *d = S[--(*Spos)];
  *c = S[--(*Spos)];
  *b = S[--(*Spos)];
  *a = S[--(*Spos)];
}

#define MIN(a, b) (a) < (b) ? (a) : (b)

#define MAX(a, b) (a) > (b) ? (a) : (b)

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
  const int16 stack_size = 4 * X.N();
  std::unique_ptr<int16[]> S = std::make_unique<int16[]>(stack_size);
  int16 Spos = 0;
  
  for (int16 i = 1; i != trees.length(); i++) {
    Spos = 0; // Empty the stack S
    
    X.CLEAR();
    ClusterTable Ti(List(trees(i)));
    Ti.TRESET();
    Ti.NVERTEX(&v, &w);
    
    do {
      if (Ti.is_leaf(&v)) {
        push(X.ENCODE(v), X.ENCODE(v), 1, 1, S, &Spos);
      } else {
        L = INF;
        R = 0;
        N = 0;
        W = 1;
        do {
          pop(&L_i, &R_i, &N_i, &W_i, S, &Spos);
          L = MIN(L, L_i);
          R = MAX(R, R_i);
          N = N + N_i;
          W = W + W_i;
          w = w - W_i;
        } while (w);
        push(L, R, N, W, S, &Spos);
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

// COMCLUSTER computes a strict consensus tree in O(knn).
// COMCLUST requires O(kn).
// trees is a list of objects of class phylo, all with the same tip labels
// (try RenumberTips(trees, trees[[1]]))
// [[Rcpp::export]]
double consensus_info (const List trees, const LogicalVector phylo) {
  
  int16 v = 0, w = 0,
    L, R, N, W,
    L_j, R_j, N_j, W_j
  ;
  
  std::vector<ClusterTable> tables;
  tables.reserve(trees.length());
  for (int16 i = trees.length(); i--; ) {
    tables.emplace_back(ClusterTable(List(trees(i))));
  }
  
  const int16
    n_trees = trees.length(),
    n_tip = tables[0].N(),
    stack_size = 4 * n_tip,
    thresh = (n_trees / 2) + 1
  ;
  
  const bool phylo_info = phylo[0];
  
  IntegerVector
    S(stack_size),
    split_count(n_tip),
    split_size(n_tip)
  ;
  
  int16
    Spos = 0,
    splits_found = 0
  ;
  
  double info = 0;
  
  // All clades in 50% consensus must occur in first 50% of trees.
  for (int16 i = 0; i != thresh; i++) {
    bool unseen_splits = false;
    for (int16 j = 1; j != n_tip - 3 + 1; j++) {
      if (!tables[i].GETSWX(&j)) {
        unseen_splits = true;
        break;
      }
    }
    if (!unseen_splits) {
      continue;
    }
    
    for (int16 j = n_tip; j--; ) {
      split_size[j] = 0;
      split_count[j] = 1; // It's in this tree!
    }
    
    for (int16 j = i + 1; j != trees.length(); j++) {
      Spos = 0; // Empty the stack S
      
      tables[i].CLEAR();
      tables[j].TRESET();
      tables[j].NVERTEX_short(&v, &w);
      int16 j_pos = 0;
      
      do {
        if (tables[j].is_leaf(&v)) {
          push(tables[i].ENCODE(v), tables[i].ENCODE(v), 1, 1, S, &Spos);
        } else {
          L = INF;
          R = 0;
          N = 0;
          W = 1;
          do {
            pop(&L_j, &R_j, &N_j, &W_j, S, &Spos);
            L = MIN(L, L_j);
            R = MAX(R, R_j);
            N = N + N_j;
            W = W + W_j;
            w = w - W_j;
          } while (w);
          push(L, R, N, W, S, &Spos);
          
          ++j_pos;
          if (tables[j].GETSWX(&j_pos)) {
            // Split has already been counted; next!
          } else {
            if (N == R - L + 1) { // L..R is contiguous, and must be tested
              if (tables[i].CLUSTONL(&L, &R)) {
                tables[j].SETSWX(&j_pos);
                assert(L > 0);
                split_count[L - 1]++;
                if (!split_size[L - 1]) split_size[L - 1] = R - L + 1;
              } else if (tables[i].CLUSTONR(&L, &R)) {
                tables[j].SETSWX(&j_pos);
                assert(R > 0);
                split_count[R - 1]++;
                if (!split_size[R - 1]) split_size[R - 1] = R - L + 1;
              }
            }
          }
        }
        tables[j].NVERTEX_short(&v, &w);
      } while (v);
    }
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
          return info;
        }
      }
    }
  }
  
  return info;
}

// [[Rcpp::export]]
IntegerVector robinson_foulds_all_pairs(List tables) {
  int16 
    v = 0, w = 0,
    L, R, N, W,
    L_i, R_i, N_i, W_i,
    n_shared,
    Spos
  ;
  const int16 n_trees = tables.length();
  if (n_trees < 2) return IntegerVector(0);
  
  IntegerVector shared(n_trees * (n_trees - 1) / 2);
  IntegerVector::iterator write_pos = shared.begin();
  
  for (int16 i = 0; i != n_trees - 1; i++) {
    Rcpp::XPtr<ClusterTable> table_i = tables(i);
    Rcpp::XPtr<ClusterTable> Xi(table_i);
    const int16 stack_size = 4 * Xi->N();
    
    for (int16 j = i + 1; j != n_trees; j++) {
      Rcpp::XPtr<ClusterTable> table_j = tables(j);
      Rcpp::XPtr<ClusterTable> Tj(table_j);
      IntegerVector S(stack_size);
      Spos = 0; // Empty the stack S
      n_shared = 0;
      
      Tj->TRESET();
      Tj->NVERTEX_short(&v, &w);
      
      do {
        if (Tj->is_leaf(&v)) {
          push(Xi->ENCODE(v), Xi->ENCODE(v), 1, 1, S, &Spos);
        } else {
          L = INF; R = 0; N = 0; W = 1;
          do {
            pop(&L_i, &R_i, &N_i, &W_i, S, &Spos);
            L = MIN(&L, &L_i);
            R = MAX(&R, &R_i);
            N = N + N_i;
            W = W + W_i;
            w = w - W_i;
          } while (w);
          push(L, R, N, W, S, &Spos);
          if (N == R - L + 1) { // L..R is contiguous, and must be tested
            if (Xi->ISCLUST(&L, &R)) n_shared++;
          }
        }
        Tj->NVERTEX_short(&v, &w); // Doesn't count all-ingroup or all-tips
      } while (v);
      *write_pos++ = n_shared;
    }
  }
  return shared;
}
