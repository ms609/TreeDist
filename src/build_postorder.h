#include <stdlib.h>

// Returns the index of the last node that it added
void add_node(int *next_node, const int node_parent, const int tips_remaining, 
              int *parent_of, int *left, int *right, const int *next_tip, const int *n_tip) {
  const int this_node = *next_node;
  parent_of[this_node] = node_parent;
  const int left_tips  = (tips_remaining > 2L) ? (rand() % (tips_remaining - 2L)) + 1L : 1L;
  const int right_tips = tips_remaining - left_tips;
  // Go right first
  if (right_tips > 1L) { // Add a node
    right[this_node - *n_tip] = ++(*next_node); // Increment value, not pointer
    add_node(next_node, this_node, right_tips, parent_of, left, right, next_tip, n_tip);
  } else { // Add a tip
    right[this_node - *n_tip] = *next_tip;
    parent_of[*next_tip] = this_node;
    next_tip++; // Increment pointer, not value. Can move this into the preceding line, [*next_tip++].
  }
  if (left_tips > 1L) {
    left[this_node - *n_tip] = ++(*next_node);
    add_node(next_node, this_node, left_tips, parent_of, left, right, next_tip, n_tip);    
  } else {
    left[this_node - *n_tip] = *next_tip;
    parent_of[*next_tip] = this_node;
    next_tip++;
  }
}

static R_NativePrimitiveArgType build_postorder_tree_t[] = {
  INTSXP, INTSXP, INTSXP, INTSXP, INTSXP
};
extern void build_postorder_tree(int *parent_of, int *left, int *right, const int *next_tip, const int *n_tip)
{
  int next_node = *n_tip;
  add_node(&next_node, *n_tip, *n_tip, parent_of, left, right, next_tip, n_tip);
}

extern SEXP BUILD_POSTORDER(SEXP ntip, SEXP tipnames) {
  // tipnames run from 0 to nTip - 1, in random order
  const int n_tip = INTEGER(ntip)[0];
  int i, *tip_names = INTEGER(tipnames);
  SEXP RESULT, PARENTOF, LEFTCHILD, RIGHTCHILD;
  PROTECT(RESULT     = allocVector(VECSXP, 3L));
  PROTECT(PARENTOF   = allocVector(INTSXP, n_tip + n_tip - 1L));
  PROTECT(LEFTCHILD  = allocVector(INTSXP, n_tip - 1L));
  PROTECT(RIGHTCHILD = allocVector(INTSXP, n_tip - 1L));
  for (i = 0; i < (n_tip + n_tip - 1L); i++) {
    INTEGER(PARENTOF)[i] = 0;
  } 
  for (i = 0; i < (n_tip - 1L); i++) {
    INTEGER(LEFTCHILD )[i] = 0;
    INTEGER(RIGHTCHILD)[i] = 0;
  }
  
  build_postorder_tree(INTEGER(PARENTOF), INTEGER(LEFTCHILD), INTEGER(RIGHTCHILD), tip_names, &n_tip);

  SET_VECTOR_ELT(RESULT, 0, PARENTOF);
  SET_VECTOR_ELT(RESULT, 1, LEFTCHILD);
  SET_VECTOR_ELT(RESULT, 2, RIGHTCHILD);
  UNPROTECT(4);
  return(RESULT);
}
