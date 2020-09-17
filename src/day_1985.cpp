#include <Rcpp.h>
using namespace Rcpp;
#include <TreeTools.h> /* for postorder_edges() */
#include <cmath> /* for log2() */
#include <memory> /* for unique_ptr */
#include "tree_distances.h"
#include "SplitList.h"

const int INF = -1;

class T {
  
  int end_of_T = 0, invocation = 0;
  static int n_edge, n_internal, n_leaves;
  
  public:
    T(IntegerMatrix);
    
    bool is_leaf(int *v) {
      return *v <= n_leaves;
    }
    
    static int edges() {
      return n_edge;
    }
    
    static int leaves() {
      return n_leaves;
    }
    
    void PREPARE(edge) {
      
    }
    
    int TEND() {
      return end_of_T;
    };
    
    void ENTER(int v, int w) {
      T(end_of_T, 0) = v;
      T(end_of_T, 1) = w;
      end_of_T++;
    }
    
    static int N() {
      return n_leaves;
    }
    
    static int M() {
      return n_internal;
    }
    
    void TRESET() {
      // This procedure prepares T for an enumeration of its entries, 
      // beginning with the first entry. 
    }
    
    void NVERTEX(int *v, int *w) {
      if (invocation != M() + N()) {
        *v = T(invocation, 0);
        *w = T(invocation, 1);
        invocation++;
      } else {
        *v = 0;
        *w = 0;
      }
    }
    
    int LEFTLEAF() {
      // If NVERTEX has returned entry <vj, wj> in T, the leftmost leaf in the
      // subtree rooted at vj has entry <vk, wk> where k = j - wj.
      // This function procedure returns Vk as its value.
      return 999;
    }
    
}

T::T(IntegerMatrix edge) {
  end_of_T = 0;
  n_edge = edge.nrow();
  IntegerVector T(n_edge);
  
}

class ClusterTable {
private:
  int enumeration = 0;
  std::unique_ptr<IntegerMatrix> X;
  
public:
  
  int ENCODE(int* v) {
    // This function procedure returns as its value the internal label 
    // assigned to leaf v
  }
  
  bool ISCLUST(int* L, int* R) {
    // This function procedure returns value true if cluster <L,R> is in X;
    // otherwise it returns value false
  }  
  
  void CLEAR() {
    // Each cluster in X has an associated switch that is either cleared or 
    // set. 
    // This procedure clears every cluster switch in X. 
  }
  
  void SETSW(int* L, int* R) {
    // If <L,R> is a cluster in X, this procedure sets the cluster switch for <L,R>. 
  }
  
  void UPDATE(){
    // This procadure inspects every cluster switch in X.
    // If the switch for cluster <L,R> is cleared, UPDATE deletes <L,R> 
    // from X; thereafter ISCLUST(X,L,R) will return the value false. 
    
  }
  
  void XRESET() {
    // This procedure prepares X for an enumeration of its clusters
    enumeration = 0;
  }
  
  void NCLUS(L, R) {
    // This procedure returns the next cluster <L,R> in the current 
    // enumeration of clusters in X.
    // If m clusters are in X, they are returned by the first m invocations 
    // of NCLUS after initialization by XRESET; thereafter NCLUS returns the
    // invalid cluster <0,0>. 
    L = X(enumeration, 0);
    R = X(enumeration, 1);
    enumeration++;
  }
}

ClusterTable::ClusterTable(T postorder_sequence) {
  
  //void BUILD(ClusterTable X) {
    // This procedure constructs in X descriptions of the clusters in a
    //  rooted tree described by the postorder sequence T with weights,
    //  BUILD assigns each leaf an internal label so that every cluster 
    //  is a set {i : L ~ i ~ R] of internal labels; thus each cluster is
    //   simply described by a pair <L,R> of internal labels, 
  X = IntegerMatrix(postorder_sequence.edges, 4);
  postorder_sequence.TRESET();
  for (int i = 1; i != N(); i++) {
    X(i, 0) = 0;
    X(i, 1) = 0;
  }
  int leafcode = 0, v, w, L, R, loc;
  
  postorder_sequence.NVERTEX(&v, &w);
  while (v) {
    if (is_leaf(&v)) {
      leafcode++;
      X(v, 2) = leafcode;
      R = leafcode;
      postorder_sequence.NVERTEX(&v, &w);
    } else {
      L = X(postorder_sequence.LEFTLEAF(), 2);
      postorder_sequence.NVERTEX(&v, &w);
      loc = w == 0 ? R : L;
      X(loc, 0) = L;
      X(loc, 1) = R;
    }
  }
  
  
}

void push (int a, int b, int c, int d, int* S) {
  *(S++) = a;
  *(S++) = b;
  *(S++) = c;
  *(S++) = d;
}

int pop (int *a, int *b, int *c, int *d, int* S) {
  *d = *(S--);
  *c = *(S--);
  *b = *(S--);
  *a = *(S--);
}

int min (int *a, int *b) {
  if (*a == INF) return *b;
  return (*a < *b ? *a : *b);
}

int max (int *a, int *b) {
  if (*a == INF) return *a;
  return (*a > *b ? *a : *b);
}

IntegerVector COMCLUST (List trees) {
  
  int *v, *w,
    *L, *R, *N, *W,
    *L_, *R_, *N_, *W_,
    *stack_0 = new int [4 * trees.length()],
    *S = stack_0;
    
  ClusterTable X = ClusterTable(trees(0)); // BUILD(T1, X)
  
  for (int i = 1; i != trees.length(); i++) {
    S = stack_0; // Empty the stack S
    
    X.CLEAR();
    Ti = T(trees(i));
    Ti.TRESET();
    Ti.NVERTEX(v, w);
    
    do {
      if (Ti.is_leaf(v)) {
        push(ENCODE(v), ENCODE(v), 1, 1, S);
      } else {
        *L = INF;
        *R = 0;
        *N = 0;
        *W = 1;
        do {
          pop(L_, R_, N_, W_, S);
          *L = min(L, L_);
          *R = max(R, R_);
          *N = *N + *N_;
          *W = *W + *W_;
          *w = *w - *W_;
        } while (*w);
        push(*L, *R, *N, *W, S);
        if (*N = *R - *L + 1) {
          SETSW(L, R);
        }
      }
      Ti.NVERTEX(v, w);
    } while (*v);
    UPDATE();
  }
  
  delete[] stack_0;
  
}
