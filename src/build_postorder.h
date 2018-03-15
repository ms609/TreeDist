#include <stdlib.h>
#include <stdio.h>
#include "RMorphy.h"

// Random number generator from http://www.cse.yorku.ca/~oz/marsaglia-rng.html
// 1+random_int%10 generates an integer from 1 to 10 [MWC renamed to random_int]
#define znew (z=36969*(z&65535)+(z>>16))
#define wnew (w=18000*(w&65535)+(w>>16))
#define random_int ((znew<<16)+wnew)

/* Global static variables: */
static unsigned long z=362436069, w=521288629;
/* Use random seeds to reset z and w*/

void insert_tip_below (const int *new_tip,
                       const int *add_below, const int *new_node,
                       int *parent_of, int *left, int *right
                       ) {
  const int old_parent = parent_of[*add_below];
  if (left[old_parent] == *add_below) {
    left[old_parent] = *new_node;
  } else {
    // The same, but on the right
    right[old_parent] = *new_node;
  }
  parent_of[*new_node] = old_parent;
  
  left[*new_node] = *new_tip;
  parent_of[*new_tip] = *new_node;
  
  right[*new_node] = *add_below;
  parent_of[*add_below] = *new_node;
}

// parent_of, left and right have been initialized with a two-taxon tree with tips 0 & 1
// left and right point n_tip _before_ left and right, so we don't need to subtract n_tip each time
// We arbitrarily choose to root our tree on tip 0, so never add to that edge or the 
// "dummy" root edge.
void build_tree(int *parent_of, int *left, int *right, const int *n_tip) {
  int tip_to_add, add_below, new_node;
  for (tip_to_add = 3; tip_to_add < *n_tip; tip_to_add++) {
    new_node = tip_to_add + *n_tip - 1;
    add_below = 1 + random_int % (tip_to_add + tip_to_add - 3); // +1 to avoid edge 0
    if (add_below < tip_to_add) { // Adding below a tip
      insert_tip_below(&tip_to_add, &add_below, &new_node, parent_of, left, right);
    } else { // Adding below an existing node
      add_below += *n_tip - tip_to_add + 1; // +1 so we never touch dummy root edge
      insert_tip_below(&tip_to_add, &add_below, &new_node, parent_of, left, right);
    }
  }
}

void move_to_node(const int *old_node_id, int *new_parent, int *new_left, int *new_right,
                  const int *old_parent, const int *old_left, const int *old_right, 
                  int *next_label, const int *n_tip) {
  const int new_node_id = *next_label;
  if (old_right[*old_node_id] > *n_tip) {
    new_right[new_node_id] = ++(*next_label);
    new_parent[*next_label] = new_node_id;
    move_to_node(&old_right[*old_node_id], new_parent, new_left, new_right,
                                 old_parent, old_left, old_right, 
                                 next_label, n_tip);    
  } else if (new_node_id != *old_node_id) { // Otherwise no change
    new_parent[old_right[*old_node_id]] = new_node_id;
    new_right[new_node_id] = old_right[*old_node_id];
  }
  if (old_left[*old_node_id] > *n_tip) {
    new_left[new_node_id] = ++(*next_label);
    new_parent[*next_label] = new_node_id;
    move_to_node(&old_left[*old_node_id], new_parent, new_left, new_right,
                                 old_parent, old_left, old_right, 
                                 next_label, n_tip);    
  } else if (new_node_id != *old_node_id) { // Otherwise no change
    new_parent[old_left[*old_node_id]] = new_node_id;
    new_left[new_node_id] = old_left[*old_node_id];
  }
}

void renumber_postorder(int *parent_of, int *left, int *right, const int *n_tip) {
  int *old_parent = malloc((*n_tip + *n_tip - 1) * sizeof(int)),
      *left_array = malloc((*n_tip - 1)          * sizeof(int)),
     *right_array = malloc((*n_tip - 1)          * sizeof(int)),
        *old_left = left_array  - *n_tip,
       *old_right = right_array - *n_tip,
       next_label = *n_tip,
                i;
  for (i = 0; i < *n_tip; i++) {
    old_parent[i] = parent_of[i];
  }
  for (i = *n_tip; i < (*n_tip + *n_tip - 1); i++) {
    old_parent[i] = parent_of[i];
    old_left  [i] = left[i];
    old_right [i] = right[i];
  }
  
  move_to_node(n_tip, parent_of, left, right, 
               old_parent, old_left, old_right, &next_label, n_tip);
  
  free(right_array);
  free(left_array);
  free(old_parent);
}

void random_tree(int *parent_of, int *left, int *right, const int *n_tip) {
  if (*n_tip < 3) {
        // Initialize with 2-tip tree
       parent_of[0] = *n_tip;
       parent_of[1] = *n_tip;
  parent_of[*n_tip] = *n_tip; // Root is its own parent
            left[0] = 0;
           right[0] = 1;
  } else {
    // Initialize with 3-tip tree, arbitrarily rooted on tip 0
             parent_of[0] = *n_tip;
             parent_of[1] = *n_tip + 1;
             parent_of[2] = *n_tip + 1;
        parent_of[*n_tip] = *n_tip; // Root is its own parent
    parent_of[*n_tip + 1] = *n_tip;
                  left[0] = 0;
                  left[1] = 1;
                 right[0] = *n_tip + 1;
                 right[1] = 2;
  }
  if (*n_tip > 3) {    
    build_tree(parent_of, left - *n_tip, right - *n_tip, n_tip);
    renumber_postorder(parent_of, left - *n_tip, right - *n_tip, n_tip);
  }
}


extern SEXP RANDOM_TREE(SEXP ntip) {
  const int n_tip = INTEGER(ntip)[0];
  SEXP RESULT = PROTECT(allocVector(VECSXP, 3)),
    PARENT_OF = PROTECT(allocVector(INTSXP, n_tip + n_tip - 1)),
         LEFT = PROTECT(allocVector(INTSXP, n_tip - 1)),
        RIGHT = PROTECT(allocVector(INTSXP, n_tip - 1));
  
  int *parent_of = INTEGER(PARENT_OF),
          *right = INTEGER(RIGHT),
           *left = INTEGER(LEFT);
 
  random_tree(parent_of, left, right, &n_tip);
  
  SET_VECTOR_ELT(RESULT, 0, PARENT_OF);
  SET_VECTOR_ELT(RESULT, 1, LEFT);
  SET_VECTOR_ELT(RESULT, 2, RIGHT);
  UNPROTECT(4);
  return(RESULT);
}

extern SEXP RANDOM_TREE_SCORE(SEXP ntip, SEXP MorphyHandl) {
  const int n_tip = INTEGER(ntip)[0];
  Morphy handl = R_ExternalPtrAddr(MorphyHandl);
  SEXP RESULT = PROTECT(allocVector(INTSXP, 1));
  int *score;
  score = INTEGER(RESULT);
  *score = 0;
  if (n_tip < 2) {
    INTEGER(RESULT)[0] = 0;
    UNPROTECT(1);
    return(RESULT);
  }
  
  // NOTE: malloc here causes segfault.
  //       Does this mean that we're accessing uninitialzied values somewhere?
  int *parent_of = calloc(n_tip + n_tip - 1 , sizeof(int)),
           *left = calloc(n_tip - 1         , sizeof(int)),
          *right = calloc(n_tip - 1         , sizeof(int));
  
  random_tree(parent_of, left, right, &n_tip);
  morphy_length(parent_of, left, right, handl, score); 
  
  free(parent_of);
  free(left);
  free(right);
  UNPROTECT(1);
  return(RESULT);
}
