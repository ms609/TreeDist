#define USE_RINTERNALS

#include <Rmath.h>
#include <math.h>
#include <R.h> 
#include <Rinternals.h>

#include "phangorn_fitch.h"
#include "reorder.h"

static const R_CMethodDef cMethods[] = {
  {"order_edges_number_nodes", (DL_FUNC) &order_edges_number_nodes, 4, order_edges_number_nodes_t},
  {"ape_neworder_phylo", (DL_FUNC) &ape_neworder_phylo, 6, ape_neworder_phylo_t},
  {"ape_neworder_pruningwise", (DL_FUNC) &ape_neworder_pruningwise, 6, ape_neworder_pruningwise_t},
  {NULL, NULL, 0, NULL}
};

static const R_CallMethodDef callMethods[] = {
  {"FITCH", (DL_FUNC) &FITCH, 8},
  {NULL, NULL, 0}
};

void R_init_myLib(DllInfo *info) {
  R_registerRoutines(info, cMethods, callMethods, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
  R_forceSymbols(info, TRUE);
}

