#include <Rcpp/Lightest>
#include <TreeTools/assert.h> // for ASSERT
#include <TreeTools/types.h> // for intx
#include <TreeTools/renumber_tree.h> // for postorder_order
using namespace Rcpp;

#define R_TO_C 1
#define IS_TIP(i) i <= n_tip

inline void add_child(const int *parent, const int *child,
                      std::unique_ptr<intx[]> & a_child,
                      std::unique_ptr<intx[]> & sibling,
                      std::unique_ptr<intx[]> & senior) {
  senior[*child] = *parent;
  // Rcout << " Parent of " << *child << " = " << (*parent) <<"\n";
  const intx existing_child = a_child[*parent];
  if (existing_child) {
    sibling[existing_child] = *child;
    sibling[*child] = existing_child;
    // Rcout << "  Siblings: " << existing_child << ", " << *child <<".\n"; 
  } else {
    a_child[*parent] = *child;
    // Rcout << "  A child of " << *parent << " = " << *child << ".\n";
  }
}

#define TODO_DELETE_RMSIB_DEBUG \
  Rcout << "sibling[" << tip << "] = sibling[parents[(tip)]] = " \
        << sibling[parents[(tip)]] << std::endl; \
  Rcout << "sibling[" << sibling[(tip)] << "] = " << (tip) << std::endl; \
  Rcout << "parents[" << tip << "] = " << parents[parents[(tip)]] << std::endl; \
  Rcout << "a_child[" << parents[(tip)] << "] = " << (tip) << std::endl;
  
// TODO might be better to re-christen as "REMOVE_TIP"
#define REMOVE_SIBLING(tip, a_child, sibling, parents)             \
  sibling[sibling[(tip)]] = 0;                                 \
  sibling[(tip)] = sibling[parents[(tip)]];                        \
  sibling[sibling[(tip)]] = (tip);                                 \
  parents[(tip)] = parents[parents[(tip)]];                        \
  a_child[parents[(tip)]] = (tip)

#define LIFT_ROOT(tip, a_child, sibling, parents)                  \
  sibling[1] = sibling[(tip)];                                     \
  parents[sibling[(tip)]] = parents[(tip)]

#define ADD_EDGE(parent, child)                                    \
  ret(*next_edge, 0) = parent;                                     \
  if (child > *n_tip) {                                            \
    ret(*next_edge, 1) = *next_node;                               \
    ++(*next_edge);                                                \
    rebuild_tree(child, next_edge, next_node, n_tip, a_child,      \
                 sibling, senior, ret);                            \
  } else {                                                         \
    ret(*next_edge, 1) = -child;                                   \
    ++(*next_edge);                                                \
  }

inline void rebuild_tree(
    const intx node,
    intx * next_edge,
    intx * next_node,
    const intx * n_tip,
    std::unique_ptr<intx[]> & a_child,
    std::unique_ptr<intx[]> & sibling,
    std::unique_ptr<intx[]> & senior,
    IntegerMatrix & ret) 
{
  const intx this_node = *next_node;
  ++(*next_node);
  
  const intx left = a_child[node];
  ADD_EDGE(this_node, left);
  ADD_EDGE(this_node, sibling[left]);
}

// edge1 and edge2 are edge matrices of binary trees with identical leaf
// labels, rooted on leaf 1, in some form of postorder.
// [[Rcpp::export]]
Rcpp::List reduce_trees(const IntegerMatrix x,
                         const IntegerMatrix y) {
  const intx
    n_edge = x.nrow(),
    n_node = n_edge / 2,
    n_tip = n_node + 1,
    n_vert = n_node + n_tip,
    root_node = n_tip + 1,
    ledger_size = n_vert + R_TO_C
  ;
  auto
    x_child_1 = std::make_unique<intx[]>(ledger_size),
    y_child_1 = std::make_unique<intx[]>(ledger_size),
    x_sibling = std::make_unique<intx[]>(ledger_size),
    y_sibling = std::make_unique<intx[]>(ledger_size),
    x_parents = std::make_unique<intx[]>(ledger_size),
    y_parents = std::make_unique<intx[]>(ledger_size)
  ;
  for (intx i = n_edge; i--; ) {
    add_child(&x(i, 0), &x(i, 1), x_child_1, x_sibling, x_parents);
    add_child(&y(i, 0), &y(i, 1), y_child_1, y_sibling, y_parents);
  }
  
  intx dropped = 0;
  for (intx it = n_tip; it--; ) {
    const intx i = it + R_TO_C;
    const intx sibling = x_sibling[i];
    // Rcout << " Cherry: " << i << ", " << sibling 
    //       << " (y = " << y_sibling[i] <<").\n";
    if (sibling && sibling <= n_tip && sibling == y_sibling[i]) {
      REMOVE_SIBLING(i, x_child_1, x_sibling, x_parents);
      REMOVE_SIBLING(i, y_child_1, y_sibling, y_parents);
      ++dropped;
      // Rcout << "   - Dropping. " << dropped << "\n";
      ++it;
    }
  }
  // Rcout << "\n dropped: " << dropped << "; ntip = " << n_tip << "\n";
  if (dropped > n_tip - 4) {
    // There's only one three-leaf topology
    return Rcpp::List::create(R_NilValue, R_NilValue);
  }
  
  // We rooted trees on leaf 1
  ASSERT(x_child_1[root_node] == 1);
  ASSERT(y_child_1[root_node] == 1);
  
  do {
    const intx x_sib = x_sibling[1];
    ASSERT(!(IS_TIP(x_sib))); // would've been collapsed
    if (IS_TIP(x_sib)) {
      // Only two leaves in tree. Not sure we can ever get here?
      // TODO delete
      return Rcpp::List::create(R_NilValue, R_NilValue);
    }
    intx x_1 = x_child_1[x_sib];
    if (!(IS_TIP(x_1))) {
      x_1 = x_sibling[x_1];
    }
    if (IS_TIP(x_1)) {
      const intx
        y_sib = y_sibling[1],
        y_1 = y_child_1[y_sib]
      ;
      if (x_1 == y_1 || x_1 == y_sibling[y_1]) {
        LIFT_ROOT(x_1, x_child_1, x_sibling, x_parents);
        LIFT_ROOT(x_1, y_child_1, y_sibling, y_parents);
        ++dropped;
      } else {
        break;
      }
    } else {
      break;
    }
  } while (true);
  
  intx
    kept_tips = n_tip - dropped,
    kept_edges = kept_tips + kept_tips - 2,
    next_edge = 0,
    next_node = kept_tips + 1
  ;
  IntegerMatrix
    x_final(kept_edges, 2),
    y_final(kept_edges, 2)
  ;
  rebuild_tree(root_node, &next_edge, &next_node, &n_tip,
               x_child_1, x_sibling, x_parents, x_final);
  next_edge = 0;
  next_node = kept_tips + 1;
  rebuild_tree(root_node, &next_edge, &next_node, &n_tip,
               y_child_1, y_sibling, y_parents, y_final);
  return Rcpp::List::create(x_final, y_final);
}
