#include <Rcpp.h>
using namespace Rcpp;
#include <TreeTools.h> /* for root_on_node() */
#include <TreeTools/root_tree.h> /* for root_on_node() */
#include <cmath> /* for log2() */
#include <memory> /* for unique_ptr */
#include "tree_distances.h"
#include "SplitList.h"

const int
  INF = INTX_MAX,
  UNINIT = -999
;

class ClusterTable {
  
  const int
    L_COL = 0,
    R_COL = 1,
    SWITCH_COL = 2,
    X_COLS = 3
  ;
  int
    n_edge,
    n_internal,
    n_leaves,
    
    n_shared = 0,
    enumeration = 0,
    v_j,
    Tlen,
    Tpos = 0
  ;
  std::unique_ptr<int[]> leftmost_leaf, T, Xarr, visit_order, decoder;
  
  public:
    ClusterTable(List); // i.e. PREPARE(T)
    
    inline bool is_leaf(const int *v) {
      return *v <= n_leaves;
    }
    
    inline const int edges() {
      return n_edge;
    }
    
    inline const int leaves() {
      return n_leaves;
    }
    
    inline int TEND() {
      return end_of_T;
    };
    
    inline void ENTER(int v, int w) {
      if (Tpos + 1 >= Tlen) std::range_error("READT T too high");
      if (Tpos < 0) std::range_error("READT T too low");
      T.get()[Tpos++] = v;
      T.get()[Tpos++] = w;
    }
    
    inline void READT(int *v, int *w) {
      *v = T[Tpos++];
      *w = T[Tpos++];
    }
    
    inline int N() {
      return n_leaves;
    }
    
    inline int M() {
      return n_internal;
    }
    
    inline void TRESET() {
      // This procedure prepares T for an enumeration of its entries, 
      // beginning with the first entry. 
      Tpos = 0;
    }
    
    inline void NVERTEX(int *v, int *w) {
      if (Tpos != Tlen) {
        if (Tpos > Tlen) throw std::range_error("Get over this programmer!");
        READT(v, w);
        v_j = *v;
      } else {
        *v = 0;
        *w = 0;
      }
    }
    
    inline int LEFTLEAF() {
      // If NVERTEX has returned entry <vj, wj> in T, the leftmost leaf in the
      // subtree rooted at vj has entry <vk, wk> where k = j - wj.
      // This function procedure returns Vk as its value.
      return leftmost_leaf[v_j];
    }
    
    inline void SET_LEFTMOST(int index, int val) {
      leftmost_leaf[index] = val;
    }
    
    inline int GET_LEFTMOST(int index) {
      return leftmost_leaf[index];
    }
    
    // Procedures to manipulate cluster tables, per Table 4 of Day 1985.
    
    inline int ENCODE(const int v) {
      // This function procedure returns as its value the internal label 
      // assigned to leaf v
      return visit_order.get()[v - 1];
    }
    
    inline int DECODE(const int v) {
      return decoder[v - 1];
    }
    
    inline int X(int row, int col) {
      return Xarr[row * X_COLS + col];
    }
    
    inline void setX(int row, int col, int value) {
      Xarr[row * X_COLS + col] = value;
    }
    
    inline bool CLUSTONL(int* L, int* R) {
      return X(*L, L_COL) == *L && X(*L, R_COL) == *R;
    }
    
    inline bool CLUSTONR(int* L, int* R) {
      return X(*R, L_COL) == *L && X(*R, R_COL) == *R;
    }
    
    inline bool ISCLUST(int* L, int* R) {
      // This function procedure returns value true if cluster <L,R> is in X;
      // otherwise it returns value false
      return CLUSTONL(L, R) || CLUSTONR(L, R);
    }  
    
    inline void CLEAR() {
      // Each cluster in X has an associated switch that is either cleared or 
      // set. 
      // This procedure clears every cluster switch in X. 
      for (int i = n_edge; i--; ) {
        Xarr[i * X_COLS + SWITCH_COL] = 0;
      }
    }
    
    inline void SETSWX(int* row) {
      setX(*row, SWITCH_COL, 1);
    }
    
    inline void SETSW(int* L, int* R) {
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
      for (int i = n_edge; i--; ) {
        int ptr = (i * X_COLS);
        if (!(Xarr[ptr + SWITCH_COL])) {
          Xarr[ptr + L_COL] = -Xarr[ptr + L_COL]; // 0
          Xarr[ptr + R_COL] = -Xarr[ptr + R_COL]; // 0
        }
      }
    }
    
    inline int SHARED() {
      return n_shared;
    }
    
    inline void ADDSHARED() {
      ++n_shared;
    }
    
    inline void XRESET() {
      // This procedure prepares X for an enumeration of its clusters
      enumeration = 0;
    }
    
    inline void NCLUS(int* L, int* R) {
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

// TODO Root tree in function; for now, must be rooted externally and in Preorder.
ClusterTable::ClusterTable(List phylo) { 
  /*
  const IntegerMatrix
    edge = phylo["edge"],
    rooted_edge = TreeTools::root_on_node(edge, 1); // Returned in preorder
  ;*/
  const IntegerMatrix edge = phylo["edge"], rooted_edge = edge;
  n_internal = phylo["Nnode"]; // = M
  CharacterVector leaf_labels = phylo["tip.label"];
  n_leaves = leaf_labels.length(); // = N
  n_edge = edge.nrow();
  Tlen = M() + N() + M() + N();
  T = std::make_unique<int[]>(Tlen);
  
  leftmost_leaf = std::make_unique<int[]>(N() + M() + 1);
  visit_order = std::make_unique<int[]>(n_leaves);
  decoder = std::make_unique<int[]>(n_leaves);
  int n_visited = 0;
  std::unique_ptr<int[]> weights = std::make_unique<int[]>(N() + M() + 1);
  
  for (int i = 1; i != n_leaves + 1; i++) {
    SET_LEFTMOST(i, i);
    weights[i] = 0;
  }
  for (int i = n_leaves + 1; i != N() + M() + 1; i++) {
    SET_LEFTMOST(i, 0);
    weights[i] = 0;
  }
  for (int i = n_edge; i--; ) {
    const int
      parent_i = rooted_edge(i, 0),
      child_i = rooted_edge(i, 1)
    ;
    if (!GET_LEFTMOST(parent_i)) SET_LEFTMOST(parent_i, GET_LEFTMOST(child_i));
    if (is_leaf(&child_i)) {
      visit_order[n_visited++] = child_i;
      decoder[child_i - 1] = n_visited;
      weights[parent_i]++;
      ENTER(child_i, 0);
    } else {
      weights[parent_i] += 1 + weights[child_i];
      ENTER(child_i, weights[child_i]);
    }
  }
  ENTER(rooted_edge(0, 0), weights[rooted_edge(0, 0)]);
  
  // BUILD Cluster table 
  Xarr = std::make_unique<int[]>(4 * edges());
  
  // This procedure constructs in X descriptions of the clusters in a
  //  rooted tree described by the postorder sequence T with weights,
  //  BUILD assigns each leaf an internal label so that every cluster 
  //  is a set {i : L ~ i ~ R] of internal labels; thus each cluster is
  //   simply described by a pair <L,R> of internal labels, 
  
  TRESET();
  for (int i = 1; i != N(); i++) {
    setX(i, L_COL, 0);
    setX(i, R_COL, 0);
  }
  int leafcode = 0, v, w, L, R = UNINIT, loc;
  
  NVERTEX(&v, &w);
  while (v) {
    if (is_leaf(&v)) {
      leafcode++;
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



inline void push (int a, int b, int c, int d, std::unique_ptr<int[]> &S, int* Spos) {
  S.get()[(*Spos)++] = a;
  S.get()[(*Spos)++] = b;
  S.get()[(*Spos)++] = c;
  S.get()[(*Spos)++] = d;
}

inline void pop (int *a, int *b, int *c, int *d, std::unique_ptr<int[]> &S, int* Spos) {
  *d = S.get()[--(*Spos)];
  *c = S.get()[--(*Spos)];
  *b = S.get()[--(*Spos)];
  *a = S.get()[--(*Spos)];
}

int min_ (int *a, int *b) {
  return (*a < *b ? *a : *b);
}

int max_ (int *a, int *b) {
  return (*a > *b ? *a : *b);
}

// COMCLUSTER computes a consensus tree in O(knn).
// COMCLUST requires O(kn).
// [[Rcpp::export]]
int COMCLUST (List trees) {
  
  int v = 0, w = 0,
    L, R, N, W,
    L_i, R_i, N_i, W_i;
  
  
  List tree_0 = trees(0);
  ClusterTable X(tree_0);
  const int stack_size = 4 * (X.N() + X.M()); // TODO this is conservative; is X.N() safe?
  std::unique_ptr<int[]> S = std::make_unique<int[]>(stack_size);
  int Spos = 0;
  
  for (int i = 1; i != trees.length(); i++) {
    Spos = 0; // Empty the stack S
    
    X.CLEAR();
    List tree_i = trees(i);
    ClusterTable Ti(tree_i);
    Ti.TRESET();
    Ti.NVERTEX(&v, &w);
    
    do {
      if (Ti.is_leaf(&v)) {
        // Rcout << Spos;
        // Rcout << " < " << stack_size;
        // check_push_safe(Spos, stack_size);
        push(Ti.ENCODE(v), Ti.ENCODE(v), 1, 1, S, &Spos);
      } else {
        L = INF;
        R = 0;
        N = 0;
        W = 1;
        do {
          pop(&L_i, &R_i, &N_i, &W_i, S, &Spos);
          L = min_(&L, &L_i);
          R = max_(&R, &R_i);
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
  
  
  IntegerMatrix ret(X.N(), 2);
  for (int i = X.N(); i--; ) {
    ret(i, 0) = X.X(i + 1, 0);
    ret(i, 1) = X.X(i + 1, 1);
  }
  return ret;
  return X.SHARED() - 2; // Subtract All-tips & All-ingroup
}
