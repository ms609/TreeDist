#define USE_RINTERNALS

#include <Rmath.h>
#include <math.h>
#include <Rinternals.h>
#include <stdlib.h> /* for NULL */
#include <R_ext/Rdynload.h>

#include "ape_reorder.h"
#include "renumber_tree.h"

#include "mpl.h"
#include "RMorphyUtils.h"
#include "RMorphy.h"
#include "build_postorder.h"


extern SEXP _TreeSearch_phangorn_bipCPP(SEXP, SEXP);

static const R_CMethodDef cMethods[] = {
  {"order_edges_number_nodes", (DL_FUNC) &order_edges_number_nodes, 3, order_edges_number_nodes_t},
  {"ape_neworder_phylo",       (DL_FUNC) &ape_neworder_phylo, 6, ape_neworder_phylo_t},
  {"ape_node_depth",           (DL_FUNC) &ape_node_depth, 7, ape_node_depth_t},
  {"ape_neworder_pruningwise", (DL_FUNC) &ape_neworder_pruningwise, 6, ape_neworder_pruningwise_t},
  {NULL, NULL, 0, NULL}
};

static const R_CallMethodDef callMethods[] = {
  {"_TreeSearch_phangorn_bipCPP",   (DL_FUNC) &_TreeSearch_phangorn_bipCPP, 2},
  {"_R_wrap_mpl_new_Morphy",        (DL_FUNC) &_R_wrap_mpl_new_Morphy, 0},
  {"_R_wrap_mpl_delete_Morphy",     (DL_FUNC) &_R_wrap_mpl_delete_Morphy, 1},
  {"_R_wrap_mpl_init_Morphy",       (DL_FUNC) &_R_wrap_mpl_init_Morphy, 3},
  {"_R_wrap_mpl_get_numtaxa",       (DL_FUNC) &_R_wrap_mpl_get_numtaxa, 1},
  {"_R_wrap_mpl_get_num_charac",    (DL_FUNC) &_R_wrap_mpl_get_num_charac, 1},
  {"_R_wrap_mpl_attach_symbols",    (DL_FUNC) &_R_wrap_mpl_attach_symbols, 2},
  {"_R_wrap_mpl_get_symbols",       (DL_FUNC) &_R_wrap_mpl_get_symbols, 1},
  {"_R_wrap_mpl_attach_rawdata",    (DL_FUNC) &_R_wrap_mpl_attach_rawdata, 2},
  {"_R_wrap_mpl_delete_rawdata",    (DL_FUNC) &_R_wrap_mpl_delete_rawdata, 1},
  {"_R_wrap_mpl_set_parsim_t",      (DL_FUNC) &_R_wrap_mpl_set_parsim_t, 3},
  {"_R_wrap_mpl_set_gaphandl",      (DL_FUNC) &_R_wrap_mpl_set_gaphandl, 2},
  {"_R_wrap_mpl_set_num_internal_nodes", (DL_FUNC) &_R_wrap_mpl_set_num_internal_nodes, 2},
  {"_R_wrap_mpl_get_num_internal_nodes", (DL_FUNC) &_R_wrap_mpl_get_num_internal_nodes, 1},
  {"_R_wrap_mpl_apply_tipdata",     (DL_FUNC) &_R_wrap_mpl_apply_tipdata, 1},
  {"_R_wrap_mpl_set_charac_weight", (DL_FUNC) &_R_wrap_mpl_set_charac_weight, 3},
  {"_R_wrap_mpl_get_charac_weight", (DL_FUNC) &_R_wrap_mpl_get_charac_weight, 2},
  {"_R_wrap_mpl_first_down_recon",  (DL_FUNC) &_R_wrap_mpl_first_down_recon, 4},
  {"_R_wrap_mpl_first_up_recon",    (DL_FUNC) &_R_wrap_mpl_first_up_recon, 5},
  {"_R_wrap_mpl_second_down_recon", (DL_FUNC) &_R_wrap_mpl_second_down_recon, 4},
  {"_R_wrap_mpl_second_up_recon",   (DL_FUNC) &_R_wrap_mpl_second_up_recon, 5},
  {"_R_wrap_mpl_update_tip",        (DL_FUNC) &_R_wrap_mpl_update_tip, 3},
  {"_R_wrap_mpl_update_lower_root", (DL_FUNC) &_R_wrap_mpl_update_lower_root, 3},
  {"MORPHYLENGTH",                  (DL_FUNC) &MORPHYLENGTH, 4},
  {"RANDOM_TREE",                   (DL_FUNC) &RANDOM_TREE, 1},
  {"RANDOM_TREE_SCORE",             (DL_FUNC) &RANDOM_TREE_SCORE, 2},

  {"RENUMBER_TREE",  (DL_FUNC) &RENUMBER_TREE,  3},
  {"RENUMBER_EDGES", (DL_FUNC) &RENUMBER_EDGES, 3},
  {NULL, NULL, 0}
};

void R_init_TreeSearch(DllInfo *dll) {
  R_registerRoutines(dll, cMethods, callMethods, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
