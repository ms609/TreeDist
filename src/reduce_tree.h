#ifndef _TREEDIST_REDUCE_TREE_H
#define _TREEDIST_REDUCE_TREE_H

#include <memory> // for unique_ptr
#include <Rcpp/Lightest>
#include <TreeTools/assert.h> // for ASSERT
#include <TreeTools/types.h> // for intx
#include <TreeTools/renumber_tree.h> // for postorder_order
using namespace Rcpp;

#define R_TO_C 1
#define IS_TIP(i) ((i) <= n_tip)
#define X_AUNT(i) x_sibling[x_parents[(i)]]
#define Y_AUNT(i) y_sibling[y_parents[(i)]]
#define Y_NIECE1(i) y_child_1[y_sibling[(i)]]
#define Y_NIECE2(i) y_sibling[Y_NIECE1(i)]
#define SAME_AUNT(x_aunt, y_aunt) x_aunt && IS_TIP(x_aunt) && x_aunt == y_aunt

extern inline void add_child(
    const int *parent, const int *child,
    std::unique_ptr<intx[]> & a_child,
    std::unique_ptr<intx[]> & sibling,
    std::unique_ptr<intx[]> & senior
  );

extern inline void rebuild_tree(
    const intx node,
    intx * next_edge,
    intx * next_node,
    const intx * n_tip,
    std::unique_ptr<intx[]> & new_no,
    std::unique_ptr<intx[]> & a_child,
    std::unique_ptr<intx[]> & sibling,
    std::unique_ptr<intx[]> & senior,
    IntegerMatrix & ret
  );

extern Rcpp::List reduce_trees(
    const IntegerMatrix x,
    const IntegerMatrix y,
    const CharacterVector original_label
  );
  
#endif