#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <float.h>
#include <R.h>
#include <Rinternals.h>
#include "mpl.h"
#include "RMorphyUtils.h"

SEXP _R_wrap_mpl_new_Morphy(void);
SEXP _R_wrap_mpl_delete_Morphy(SEXP MorphyHandl);
SEXP _R_wrap_mpl_init_Morphy(SEXP Rntax, SEXP Rnchar, SEXP MorphHandl);    
SEXP _R_wrap_mpl_get_numtaxa(SEXP MorphHandl);
SEXP _R_wrap_mpl_get_num_charac(SEXP MorphHandl);
SEXP _R_wrap_mpl_attach_symbols(SEXP Rsymbols, SEXP MorphyHandl);
SEXP _R_wrap_mpl_get_symbols(SEXP MorphyHandl);
SEXP _R_wrap_mpl_attach_rawdata(SEXP Rmatrix, SEXP MorphyHandl);
SEXP _R_wrap_mpl_delete_rawdata(SEXP MorphyHandl);
SEXP _R_wrap_mpl_set_parsim_t(SEXP RcharID, SEXP Rchtype, SEXP MorphyHandl);
SEXP _R_wrap_mpl_set_gaphandl(SEXP Rgaptype, SEXP MorphyHandl);
SEXP _R_wrap_mpl_set_num_internal_nodes(SEXP Rnnodes, SEXP MorphyHandl);
SEXP _R_wrap_mpl_get_num_internal_nodes(SEXP MorphyHandl);
SEXP _R_wrap_mpl_apply_tipdata(SEXP MorphyHandl);
SEXP _R_wrap_mpl_set_charac_weight(SEXP RcharID, SEXP Rweight, SEXP MorphyHandl);
SEXP _R_wrap_mpl_get_charac_weight(SEXP RcharID, SEXP MorphyHandl);
SEXP _R_wrap_mpl_first_down_recon(SEXP Rnode_id, SEXP Rleft_id, SEXP Rright_id, SEXP MorphyHandl);
SEXP _R_wrap_mpl_first_up_recon(SEXP Rnode_id, SEXP Rleft_id, SEXP Rright_id, SEXP Ranc_id, SEXP MorphyHandl);
SEXP _R_wrap_mpl_second_down_recon(SEXP Rnode_id, SEXP Rleft_id, SEXP Rright_id, SEXP MorphyHandl);
SEXP _R_wrap_mpl_second_up_recon(SEXP Rnode_id, SEXP Rleft_id, SEXP Rright_id, SEXP Ranc_id, SEXP MorphyHandl);
SEXP _R_wrap_mpl_update_tip(SEXP tip_id, SEXP anc_id, SEXP MorphyHandl);
SEXP _R_wrap_mpl_update_lower_root(SEXP lower_id, SEXP upper_id, SEXP MorphyHandl);
void morphy_length(const int *ancestor, const int *left, const int *right, Morphy handl, int *score);
SEXP MORPHYLENGTH(SEXP R_ancestors, SEXP R_left, SEXP R_right, SEXP MorphyHandl);