#include <Rcpp/Lightest>
#include <TreeTools/assert.h> // for ASSERT
#include <TreeTools/types.h> // for intx
#include <TreeTools/renumber_tree.h> // for postorder_order
using namespace Rcpp;
// #define TD_DEBUG

#define R_TO_C 1
#define IS_TIP(i) ((i) <= n_tip)
#define X_AUNT(i) x_sibling[x_parents[(i)]]
#define Y_AUNT(i) y_sibling[y_parents[(i)]]
#define Y_NIECE1(i) y_child_1[y_sibling[(i)]]
#define Y_NIECE2(i) y_sibling[Y_NIECE1(i)]
#define SAME_AUNT(x_aunt, y_aunt) x_aunt && IS_TIP(x_aunt) && x_aunt == y_aunt

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
  Rcout << "sibling[" << tip << "] <- sibling[parents[(tip)]] = " \
        << sibling[parents[(tip)]] << std::endl; \
  Rcout << "sibling[" << sibling[(tip)] << "] = " << (tip) << std::endl; \
  Rcout << "parents[" << tip << "] = " << parents[parents[(tip)]] << std::endl; \
  Rcout << "a_child[" << parents[(tip)] << "] = " << (tip) << std::endl;
  
// TODO might be better to re-christen as "REMOVE_TIP"
#define REMOVE_SIBLING(tip, a_child, sibling, parents)             \
  sibling[sibling[(tip)]] = 0;                                     \
  sibling[(tip)] = sibling[parents[(tip)]];                        \
  sibling[sibling[(tip)]] = (tip);                                 \
  parents[(tip)] = parents[parents[(tip)]];                        \
  a_child[parents[(tip)]] = (tip)

#define TODO_DELETE_RMTIP_DEBUG                                \
Rcout << "\n\n ==== Remove tip " << tip << ". ====\n";       \
Rcout << "parents[" << sibling[(tip)] << "] <- parents[" << parents[(tip)] \
      << "] = " << parents[parents[(tip)]] << ";\n";         \
Rcout << "sibling[" << sibling[parents[(tip)]] << "] = " << sibling[(tip)] << ";\n";\
Rcout << "sibling["<<sibling[(tip)]<<"] = "<<sibling[parents[(tip)]]<<"\n"; \
Rcout << "a_child["<<parents[sibling[(tip)]]<<"] = "<< sibling[sibling[(tip)]] <<"\n";\
Rcout << "";
  
#define REMOVE_TIP(tip, a_child, sibling, parents)           \
  ASSERT(parents[parents[(tip)]])                            \
  parents[sibling[(tip)]] = parents[parents[(tip)]];         \
  sibling[sibling[parents[(tip)]]] = sibling[(tip)];         \
  sibling[sibling[(tip)]] = sibling[parents[(tip)]];         \
  a_child[parents[sibling[(tip)]]] = sibling[sibling[(tip)]];         \
  sibling[(tip)] = 0                                           
  
#define REDUCE_CHAIN                                           \
  ASSERT(gg_aunt > 1);                                         \
  REMOVE_TIP(gg_aunt, x_child_1, x_sibling, x_parents);        \
  REMOVE_TIP(gg_aunt, y_child_1, y_sibling, y_parents);        \
  ++dropped

#define TODO_DELETE_LIFT_ROOT_DEBUG \
  Rcout << "   ^ Lifting; sibling[1] = " << sibling[(tip)];    \
  Rcout << "; sibling[" << sibling[(tip)] << "] = 1 \n";       \
  Rcout << "     parents[" << sibling[(tip)] << "] = " << parents[1] << "\n";\
  Rcout << ""
  
#define LIFT_ROOT(tip, a_child, sibling, parents)              \
  sibling[1] = sibling[(tip)];                                 \
  sibling[sibling[(tip)]] = 1;                                 \
  parents[sibling[(tip)]] = parents[1];                        \
  ASSERT(a_child[parents[1]] == 1);                            \
  sibling[(tip)] = 0                            

#define ADD_EDGE(parent, child)                                \
  --(*next_edge);                                              \
  ASSERT(*next_edge >= 0);                                     \
  ASSERT(*next_edge < ret.nrow());                             \
  ret(*next_edge, 0) = parent;                                 \
  if (child > *n_tip) {                                        \
    ret(*next_edge, 1) = *next_node;                           \
    rebuild_tree(child, next_edge, next_node, n_tip, new_no,   \
                 a_child, sibling, senior, ret);               \
  } else {                                                     \
    ret(*next_edge, 1) = new_no[child];                        \
  }

inline void rebuild_tree(
    const intx node,
    intx * next_edge,
    intx * next_node,
    const intx * n_tip,
    std::unique_ptr<intx[]> & new_no,
    std::unique_ptr<intx[]> & a_child,
    std::unique_ptr<intx[]> & sibling,
    std::unique_ptr<intx[]> & senior,
    IntegerMatrix & ret) 
{
  const intx this_node = *next_node;
  ++(*next_node);
 
  const intx left = a_child[node];
#ifdef TD_DEBUG
  Rcout << "  r Rebuilding " << node << " -> {" << left <<  ", "
        << sibling[left] << "}.\n";
#endif
  ASSERT(left);
  ASSERT(sibling[left]);
  
  ADD_EDGE(this_node, left);
  ADD_EDGE(this_node, sibling[left]);
#ifdef TD_DEBUG
  Rcout << "  x Done with " << node << " -> {" << left <<  ", "
        << sibling[left] << "}.\n";
#endif
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
  // _sibling = 0 denotes that a leaf has been removed from a tree.
  
  for (intx i = n_edge; i--; ) {
    add_child(&x(i, 0), &x(i, 1), x_child_1, x_sibling, x_parents);
    add_child(&y(i, 0), &y(i, 1), y_child_1, y_sibling, y_parents);
  }
  
  intx dropped = 0;
  for (intx it = n_tip; it--; ) {
    const intx i = it + R_TO_C;
    const intx sibling = x_sibling[i];
    if (sibling && sibling <= n_tip && sibling == y_sibling[i]) {
      REMOVE_SIBLING(i, x_child_1, x_sibling, x_parents);
      REMOVE_SIBLING(i, y_child_1, y_sibling, y_parents);
      ++dropped;
      ++it;
    }
  }
  
#ifdef TD_DEBUG
  Rcout << "\n Collapsed " << dropped << " cherries; ntip0 = " << n_tip << "\n";
#endif
  if (dropped > n_tip - 4) {
    // There's only one three-leaf topology
    return Rcpp::List::create(R_NilValue, R_NilValue);
  }
  
  // We rooted trees on leaf 1
  ASSERT(x_child_1[root_node] == 1 || x_sibling[x_child_1[root_node]] == 1);
  ASSERT(y_child_1[root_node] == 1 || y_sibling[y_child_1[root_node]] == 1);
  
  do {
    const intx x_sib = x_sibling[1];
    ASSERT(x_sib);
    ASSERT(!(IS_TIP(x_sib))); // would've been collapsed
    intx x_1 = x_child_1[x_sib];
#ifdef TD_DEBUG
    Rcout << " ^ 1-sibling " << x_sib << ", whose children = {";
#endif
    if (!(IS_TIP(x_1))) {
#ifdef TD_DEBUG
      Rcout << x_1 << ", **" << x_sibling[x_1] << "**}.\n";
#endif
      x_1 = x_sibling[x_1];
    } else {
#ifdef TD_DEBUG
      Rcout << "**" <<  x_1 << "**, " << x_sibling[x_1] << "}.\n";
#endif
    }
    if (IS_TIP(x_1)) {
      const intx
        y_sib = y_sibling[1],
        y_1 = y_child_1[y_sib]
      ;
#ifdef TD_DEBUG
      Rcout << " ^ Checking for match for " << x_1 << " in Y: {" 
            << y_1 << ", " << y_sibling[y_1] << "}\n";
#endif
      
      if (x_1 == y_1 || x_1 == y_sibling[y_1]) {
#ifdef TD_DEBUG
        Rcout << "  ^ Lifting root by dropping " << x_1 << ".\n";
#endif
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
  
#ifdef TD_DEBUG
  Rcout << "\n ^ Can't lift root further. \n";
#endif
  
  if (dropped > n_tip - 4) {
    // There's only one three-leaf topology
    return Rcpp::List::create(R_NilValue, R_NilValue);
  }
  
#ifdef TD_DEBUG
  Rcout << "\n == Chain rule == \n";
#endif
  
  for (intx it = n_tip - 1; it--; ) {
    // 1 can't be at the head of a chain: it's the root.
    const int i = it + R_TO_C + 1;
    if (!x_sibling[i]) {
#ifdef TD_DEBUG
      Rcout << " . Tip " << i << " already expunged.\n";
#endif
      continue;
    }
    const int aunt = X_AUNT(i);
    if (!aunt || !IS_TIP(aunt)) {
#ifdef TD_DEBUG
      Rcout << "> No chain at " << i << ":  aunt = " << aunt << ".\n";
#endif
      continue;
    }
    const int g_aunt = X_AUNT(aunt);
    if (!g_aunt || !IS_TIP(g_aunt)) {
#ifdef TD_DEBUG
      Rcout << "> No chain at " << i << "-" << aunt
            << ":  gt_aunt = " << g_aunt << ".\n";
#endif
      continue;
    }
    int gg_aunt = X_AUNT(g_aunt);
#ifdef TD_DEBUG
    Rcout << " o Candidate chain: " << i << "." << aunt << "."
          << g_aunt << "-(" << gg_aunt;
#endif
    if (aunt == Y_AUNT(i) &&
        g_aunt == Y_AUNT(aunt)) {
      // Case 1: Same direction
      while(SAME_AUNT(gg_aunt, Y_AUNT(g_aunt))) {
#ifdef TD_DEBUG
        Rcout << "!x" << gg_aunt;
#endif
        REDUCE_CHAIN;
        gg_aunt = X_AUNT(g_aunt);
      }
    } else if (i == Y_AUNT(aunt) &&
      aunt == Y_AUNT(g_aunt)) {
      // Case 2: Opposite direction
#ifdef TD_DEBUG
      Rcout << " = " << Y_NIECE1(g_aunt)
            << "|" << Y_NIECE2(g_aunt)
            << "|" << y_sibling[g_aunt];
#endif
      while(gg_aunt &&
      IS_TIP(gg_aunt) && (
          gg_aunt == Y_NIECE1(g_aunt) ||
            gg_aunt == Y_NIECE2(g_aunt) ||
            gg_aunt == y_sibling[g_aunt])) {
#ifdef TD_DEBUG
        Rcout << ")- rm " << x_sibling[x_parents[g_aunt]] << " -(";
#endif
        REDUCE_CHAIN;
        gg_aunt = X_AUNT(g_aunt);
      }
    }
#ifdef TD_DEBUG
    Rcout << ").\n";
#endif
  }
  
  
  auto new_no = std::make_unique<intx[]>(ledger_size);
  const intx kept_tips = n_tip - dropped;
  intx tips_left = kept_tips;
  LogicalVector tip_kept(n_tip);
  for (intx i = n_tip; i--; ) {
    if (x_sibling[i + R_TO_C]) {
      tip_kept[i] = true;
      new_no[i + R_TO_C] = tips_left;
      --tips_left;
    }
  }
  const intx kept_edges = kept_tips + kept_tips - 2;
  intx
    next_edge = kept_edges,
    next_node = kept_tips + 1
  ;
  IntegerMatrix
    x_final(kept_edges, 2),
    y_final(kept_edges, 2)
  ;
  rebuild_tree(root_node, &next_edge, &next_node, &n_tip, new_no,
               x_child_1, x_sibling, x_parents, x_final);
#ifdef TD_DEBUG
  Rcout << "\n\n == Now to rebuild tree 2 ==\n";
#endif
  ASSERT(next_node == kept_tips + kept_tips);
  ASSERT(next_edge == 0);
  next_edge = kept_edges;
  next_node = kept_tips + 1;
  rebuild_tree(root_node, &next_edge, &next_node, &n_tip, new_no,
               y_child_1, y_sibling, y_parents, y_final);
  ASSERT(next_node == kept_tips + kept_tips);
  ASSERT(next_edge == 0);
  
  return Rcpp::List::create(x_final, y_final, tip_kept);
}
