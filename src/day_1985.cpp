#include <Rcpp.h>
using namespace Rcpp;
#include <TreeTools.h> /* for root_on_node() */
#include <TreeTools/root_tree.h> /* for root_on_node() */
#include <cmath> /* for log2() */
#include <memory> /* for unique_ptr */
#include "tree_distances.h"
#include "SplitList.h"

const int INF = INTX_MAX;

class ClusterTable {
  
  const int
    L_COL = 0,
    R_COL = 1,
    SWITCH_COL = 2,
    X_COLS = 3
  ;
  int n_edge, n_internal, n_leaves;
  int 
    end_of_T = 0,
    invocation = 0,
    enumeration = 0,
    v_j,
    *Xarr,
    *leftmost_leaf, *visit_order, *decoder;
  IntegerMatrix T;
  
  public:
    ClusterTable(List); // i.e. PREPARE(T)
    ~ClusterTable();
    
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
    
    void ENTER(int v, int w) {
      T(end_of_T, 0) = v;
      T(end_of_T, 1) = w;
      end_of_T++;
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
      invocation = 0;
    }
    
    inline void NVERTEX(int *v, int *w) {
      if (invocation != M() + N()) {
        *v = T(invocation, 0);
        *w = T(invocation, 1);
        invocation++;
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
    
    // Procedures to manipulate cluster tables, per Table 4 of Day 1985.
    
    inline int ENCODE(const int* v) {
      // This function procedure returns as its value the internal label 
      // assigned to leaf v
      return visit_order[*v - 1];
    }
    
    inline int DECODE(const int* v) {
      return decoder[*v - 1];
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
        SETSWX(L);
      } else if (CLUSTONR(L, R)) {
        SETSWX(R);
      }
    }
    
    inline void UPDATE(){
      // This procadure inspects every cluster switch in X.
      // If the switch for cluster <L,R> is cleared, UPDATE deletes <L,R> 
      // from X; thereafter ISCLUST(X,L,R) will return the value false. 
      for (int i = n_edge; i--; ) {
        int *ptr = Xarr + (i * X_COLS);
        if (*(ptr + SWITCH_COL)) {
          *(ptr + L_COL) = 0;
          *(ptr + R_COL) = 0;
        }
      }
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
  IntegerMatrix T(n_edge, 2);
  end_of_T = 0;
  
  leftmost_leaf = new int[n_leaves + 1];
  visit_order = new int[n_leaves];
  decoder = new int[n_leaves];
  int n_visited = 0;
  int* weights = new int[N() + M() + 1]; 
  
  for (int i = 1; i != n_leaves; i++) {
    leftmost_leaf[i] = i;
    weights[i] = 0;
  }
  for (int i = n_edge; i--; ) {
    const int
      parent_i = rooted_edge(i, 0),
      child_i = rooted_edge(i, 1)
    ;
    leftmost_leaf[parent_i] = leftmost_leaf[child_i];
    if (this->is_leaf(&child_i)) {
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
  delete[] weights;
  
  // BUILD Cluster table 
  Xarr = new int[4 * edges()];
  
  // This procedure constructs in X descriptions of the clusters in a
  //  rooted tree described by the postorder sequence T with weights,
  //  BUILD assigns each leaf an internal label so that every cluster 
  //  is a set {i : L ~ i ~ R] of internal labels; thus each cluster is
  //   simply described by a pair <L,R> of internal labels, 
  
  TRESET();
  for (int i = 1; i != N(); i++) {
    setX(i, 0, 0);
    setX(i, 1, 0);
  }
  int leafcode = 0, *v, *w, *L, *R, loc;
  
  NVERTEX(v, w);
  while (*v) {
    if (is_leaf(v)) {
      leafcode++;
      setX(*v, 2, leafcode);
      *R = leafcode;
      NVERTEX(v, w);
    } else {
      *L = X(LEFTLEAF(), 2);
      NVERTEX(v, w);
      loc = *w == 0 ? *R : *L;
      setX(loc, 0, *L);
      setX(loc, 1, *R);
    }
  }
}

ClusterTable::~ClusterTable() {
  delete[] leftmost_leaf;
  delete[] visit_order;
  delete[] decoder;
  delete[] Xarr;
}


void push (int a, int b, int c, int d, int* S) {
  *(S++) = a;
  *(S++) = b;
  *(S++) = c;
  *(S++) = d;
}

void pop (int *a, int *b, int *c, int *d, int* S) {
  *d = *(S--);
  *c = *(S--);
  *b = *(S--);
  *a = *(S--);
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
IntegerMatrix COMCLUST (List trees) {
  
  int *v = 0, *w = 0,
    *L, *R, *N, *W,
    *L_, *R_, *N_, *W_,
    *stack_0 = new int [4 * trees.length()],
    *S = stack_0;
  
  List tree_0 = trees(0);
  ClusterTable X(tree_0);
  
  for (int i = 1; i != trees.length(); i++) {
    S = stack_0; // Empty the stack S
    
    X.CLEAR();
    List tree_i = trees(i);
    ClusterTable Ti(tree_i);
    Ti.TRESET();
    Ti.NVERTEX(v, w);
    
    do {
      if (Ti.is_leaf(v)) {
        push(Ti.ENCODE(v), Ti.ENCODE(v), 1, 1, S);
      } else {
        *L = INF;
        *R = 0;
        *N = 0;
        *W = 1;
        do {
          pop(L_, R_, N_, W_, S);
          *L = min_(L, L_);
          *R = max_(R, R_);
          *N = *N + *N_;
          *W = *W + *W_;
          *w = *w - *W_;
        } while (*w);
        push(*L, *R, *N, *W, S);
        if (*N == *R - *L + 1) { // L..R is contiguous, and must be tested
          X.SETSW(L, R);
        }
      }
      Ti.NVERTEX(v, w);
    } while (*v);
    X.UPDATE();
  }
  
  delete[] stack_0;
  
  IntegerMatrix ret(X.N(), 2);
  for (int i = X.N(); i--; ) {
    ret(i, 0) = X.X(i + 1, 0);
    ret(i, 1) = X.X(i + 1, 1);
  }
  return ret;
}
