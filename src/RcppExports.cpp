// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// ClusterTable_new
SEXP ClusterTable_new(List phylo);
RcppExport SEXP _TreeDist_ClusterTable_new(SEXP phyloSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type phylo(phyloSEXP);
    rcpp_result_gen = Rcpp::wrap(ClusterTable_new(phylo));
    return rcpp_result_gen;
END_RCPP
}
// ClusterTable_matrix
IntegerMatrix ClusterTable_matrix(SEXP xp);
RcppExport SEXP _TreeDist_ClusterTable_matrix(SEXP xpSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type xp(xpSEXP);
    rcpp_result_gen = Rcpp::wrap(ClusterTable_matrix(xp));
    return rcpp_result_gen;
END_RCPP
}
// ClusterTable_decode
IntegerVector ClusterTable_decode(SEXP xp);
RcppExport SEXP _TreeDist_ClusterTable_decode(SEXP xpSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type xp(xpSEXP);
    rcpp_result_gen = Rcpp::wrap(ClusterTable_decode(xp));
    return rcpp_result_gen;
END_RCPP
}
// COMCLUST
int COMCLUST(List trees);
RcppExport SEXP _TreeDist_COMCLUST(SEXP treesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type trees(treesSEXP);
    rcpp_result_gen = Rcpp::wrap(COMCLUST(trees));
    return rcpp_result_gen;
END_RCPP
}
// consensus_info
double consensus_info(const List trees, const LogicalVector phylo, const NumericVector p);
RcppExport SEXP _TreeDist_consensus_info(SEXP treesSEXP, SEXP phyloSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List >::type trees(treesSEXP);
    Rcpp::traits::input_parameter< const LogicalVector >::type phylo(phyloSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(consensus_info(trees, phylo, p));
    return rcpp_result_gen;
END_RCPP
}
// robinson_foulds_all_pairs
IntegerVector robinson_foulds_all_pairs(List tables);
RcppExport SEXP _TreeDist_robinson_foulds_all_pairs(SEXP tablesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type tables(tablesSEXP);
    rcpp_result_gen = Rcpp::wrap(robinson_foulds_all_pairs(tables));
    return rcpp_result_gen;
END_RCPP
}
// lapjv
List lapjv(NumericMatrix x, NumericVector maxX);
RcppExport SEXP _TreeDist_lapjv(SEXP xSEXP, SEXP maxXSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type maxX(maxXSEXP);
    rcpp_result_gen = Rcpp::wrap(lapjv(x, maxX));
    return rcpp_result_gen;
END_RCPP
}
// cpp_mast
int cpp_mast(IntegerMatrix edge1, IntegerMatrix edge2, IntegerVector nTip);
RcppExport SEXP _TreeDist_cpp_mast(SEXP edge1SEXP, SEXP edge2SEXP, SEXP nTipSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type edge1(edge1SEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type edge2(edge2SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nTip(nTipSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_mast(edge1, edge2, nTip));
    return rcpp_result_gen;
END_RCPP
}
// cpp_nni_distance
IntegerVector cpp_nni_distance(const IntegerMatrix edge1, const IntegerMatrix edge2, const IntegerVector nTip);
RcppExport SEXP _TreeDist_cpp_nni_distance(SEXP edge1SEXP, SEXP edge2SEXP, SEXP nTipSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerMatrix >::type edge1(edge1SEXP);
    Rcpp::traits::input_parameter< const IntegerMatrix >::type edge2(edge2SEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type nTip(nTipSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_nni_distance(edge1, edge2, nTip));
    return rcpp_result_gen;
END_RCPP
}
// path_vector
IntegerVector path_vector(IntegerMatrix edge);
RcppExport SEXP _TreeDist_path_vector(SEXP edgeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type edge(edgeSEXP);
    rcpp_result_gen = Rcpp::wrap(path_vector(edge));
    return rcpp_result_gen;
END_RCPP
}
// vec_diff_euclidean
NumericMatrix vec_diff_euclidean(const IntegerMatrix vec1, const IntegerMatrix vec2);
RcppExport SEXP _TreeDist_vec_diff_euclidean(SEXP vec1SEXP, SEXP vec2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerMatrix >::type vec1(vec1SEXP);
    Rcpp::traits::input_parameter< const IntegerMatrix >::type vec2(vec2SEXP);
    rcpp_result_gen = Rcpp::wrap(vec_diff_euclidean(vec1, vec2));
    return rcpp_result_gen;
END_RCPP
}
// pair_diff_euclidean
NumericVector pair_diff_euclidean(const IntegerMatrix vecs);
RcppExport SEXP _TreeDist_pair_diff_euclidean(SEXP vecsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerMatrix >::type vecs(vecsSEXP);
    rcpp_result_gen = Rcpp::wrap(pair_diff_euclidean(vecs));
    return rcpp_result_gen;
END_RCPP
}
// path_vector2
IntegerVector path_vector2(IntegerMatrix edge);
RcppExport SEXP _TreeDist_path_vector2(SEXP edgeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type edge(edgeSEXP);
    rcpp_result_gen = Rcpp::wrap(path_vector2(edge));
    return rcpp_result_gen;
END_RCPP
}
// reduce_trees
Rcpp::List reduce_trees(const IntegerMatrix x, const IntegerMatrix y, const CharacterVector original_label);
RcppExport SEXP _TreeDist_reduce_trees(SEXP xSEXP, SEXP ySEXP, SEXP original_labelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< const IntegerMatrix >::type y(ySEXP);
    Rcpp::traits::input_parameter< const CharacterVector >::type original_label(original_labelSEXP);
    rcpp_result_gen = Rcpp::wrap(reduce_trees(x, y, original_label));
    return rcpp_result_gen;
END_RCPP
}
// mismatch_size
IntegerVector mismatch_size(const RawMatrix x, const RawMatrix y);
RcppExport SEXP _TreeDist_mismatch_size(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const RawMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< const RawMatrix >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(mismatch_size(x, y));
    return rcpp_result_gen;
END_RCPP
}
// confusion
IntegerVector confusion(const RawMatrix x, const RawMatrix y);
RcppExport SEXP _TreeDist_confusion(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const RawMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< const RawMatrix >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(confusion(x, y));
    return rcpp_result_gen;
END_RCPP
}
// keep_and_reroot
List keep_and_reroot(const List tree1, const List tree2, const LogicalVector keep);
RcppExport SEXP _TreeDist_keep_and_reroot(SEXP tree1SEXP, SEXP tree2SEXP, SEXP keepSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List >::type tree1(tree1SEXP);
    Rcpp::traits::input_parameter< const List >::type tree2(tree2SEXP);
    Rcpp::traits::input_parameter< const LogicalVector >::type keep(keepSEXP);
    rcpp_result_gen = Rcpp::wrap(keep_and_reroot(tree1, tree2, keep));
    return rcpp_result_gen;
END_RCPP
}
// keep_and_reduce
List keep_and_reduce(const List tree1, const List tree2, const LogicalVector keep);
RcppExport SEXP _TreeDist_keep_and_reduce(SEXP tree1SEXP, SEXP tree2SEXP, SEXP keepSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List >::type tree1(tree1SEXP);
    Rcpp::traits::input_parameter< const List >::type tree2(tree2SEXP);
    Rcpp::traits::input_parameter< const LogicalVector >::type keep(keepSEXP);
    rcpp_result_gen = Rcpp::wrap(keep_and_reduce(tree1, tree2, keep));
    return rcpp_result_gen;
END_RCPP
}
// cpp_robinson_foulds_distance
List cpp_robinson_foulds_distance(const RawMatrix x, const RawMatrix y, const IntegerVector nTip);
RcppExport SEXP _TreeDist_cpp_robinson_foulds_distance(SEXP xSEXP, SEXP ySEXP, SEXP nTipSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const RawMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< const RawMatrix >::type y(ySEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type nTip(nTipSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_robinson_foulds_distance(x, y, nTip));
    return rcpp_result_gen;
END_RCPP
}
// cpp_robinson_foulds_info
List cpp_robinson_foulds_info(const RawMatrix x, const RawMatrix y, const IntegerVector nTip);
RcppExport SEXP _TreeDist_cpp_robinson_foulds_info(SEXP xSEXP, SEXP ySEXP, SEXP nTipSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const RawMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< const RawMatrix >::type y(ySEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type nTip(nTipSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_robinson_foulds_info(x, y, nTip));
    return rcpp_result_gen;
END_RCPP
}
// cpp_matching_split_distance
List cpp_matching_split_distance(const RawMatrix x, const RawMatrix y, const IntegerVector nTip);
RcppExport SEXP _TreeDist_cpp_matching_split_distance(SEXP xSEXP, SEXP ySEXP, SEXP nTipSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const RawMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< const RawMatrix >::type y(ySEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type nTip(nTipSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_matching_split_distance(x, y, nTip));
    return rcpp_result_gen;
END_RCPP
}
// cpp_jaccard_similarity
List cpp_jaccard_similarity(const RawMatrix x, const RawMatrix y, const IntegerVector nTip, const NumericVector k, const LogicalVector allowConflict);
RcppExport SEXP _TreeDist_cpp_jaccard_similarity(SEXP xSEXP, SEXP ySEXP, SEXP nTipSEXP, SEXP kSEXP, SEXP allowConflictSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const RawMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< const RawMatrix >::type y(ySEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type nTip(nTipSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type k(kSEXP);
    Rcpp::traits::input_parameter< const LogicalVector >::type allowConflict(allowConflictSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_jaccard_similarity(x, y, nTip, k, allowConflict));
    return rcpp_result_gen;
END_RCPP
}
// cpp_msi_distance
List cpp_msi_distance(const RawMatrix x, const RawMatrix y, const IntegerVector nTip);
RcppExport SEXP _TreeDist_cpp_msi_distance(SEXP xSEXP, SEXP ySEXP, SEXP nTipSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const RawMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< const RawMatrix >::type y(ySEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type nTip(nTipSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_msi_distance(x, y, nTip));
    return rcpp_result_gen;
END_RCPP
}
// cpp_mutual_clustering
List cpp_mutual_clustering(const RawMatrix x, const RawMatrix y, const IntegerVector nTip);
RcppExport SEXP _TreeDist_cpp_mutual_clustering(SEXP xSEXP, SEXP ySEXP, SEXP nTipSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const RawMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< const RawMatrix >::type y(ySEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type nTip(nTipSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_mutual_clustering(x, y, nTip));
    return rcpp_result_gen;
END_RCPP
}
// cpp_shared_phylo
List cpp_shared_phylo(const RawMatrix x, const RawMatrix y, const IntegerVector nTip);
RcppExport SEXP _TreeDist_cpp_shared_phylo(SEXP xSEXP, SEXP ySEXP, SEXP nTipSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const RawMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< const RawMatrix >::type y(ySEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type nTip(nTipSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_shared_phylo(x, y, nTip));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_TreeDist_ClusterTable_new", (DL_FUNC) &_TreeDist_ClusterTable_new, 1},
    {"_TreeDist_ClusterTable_matrix", (DL_FUNC) &_TreeDist_ClusterTable_matrix, 1},
    {"_TreeDist_ClusterTable_decode", (DL_FUNC) &_TreeDist_ClusterTable_decode, 1},
    {"_TreeDist_COMCLUST", (DL_FUNC) &_TreeDist_COMCLUST, 1},
    {"_TreeDist_consensus_info", (DL_FUNC) &_TreeDist_consensus_info, 3},
    {"_TreeDist_robinson_foulds_all_pairs", (DL_FUNC) &_TreeDist_robinson_foulds_all_pairs, 1},
    {"_TreeDist_lapjv", (DL_FUNC) &_TreeDist_lapjv, 2},
    {"_TreeDist_cpp_mast", (DL_FUNC) &_TreeDist_cpp_mast, 3},
    {"_TreeDist_cpp_nni_distance", (DL_FUNC) &_TreeDist_cpp_nni_distance, 3},
    {"_TreeDist_path_vector", (DL_FUNC) &_TreeDist_path_vector, 1},
    {"_TreeDist_vec_diff_euclidean", (DL_FUNC) &_TreeDist_vec_diff_euclidean, 2},
    {"_TreeDist_pair_diff_euclidean", (DL_FUNC) &_TreeDist_pair_diff_euclidean, 1},
    {"_TreeDist_path_vector2", (DL_FUNC) &_TreeDist_path_vector2, 1},
    {"_TreeDist_reduce_trees", (DL_FUNC) &_TreeDist_reduce_trees, 3},
    {"_TreeDist_mismatch_size", (DL_FUNC) &_TreeDist_mismatch_size, 2},
    {"_TreeDist_confusion", (DL_FUNC) &_TreeDist_confusion, 2},
    {"_TreeDist_keep_and_reroot", (DL_FUNC) &_TreeDist_keep_and_reroot, 3},
    {"_TreeDist_keep_and_reduce", (DL_FUNC) &_TreeDist_keep_and_reduce, 3},
    {"_TreeDist_cpp_robinson_foulds_distance", (DL_FUNC) &_TreeDist_cpp_robinson_foulds_distance, 3},
    {"_TreeDist_cpp_robinson_foulds_info", (DL_FUNC) &_TreeDist_cpp_robinson_foulds_info, 3},
    {"_TreeDist_cpp_matching_split_distance", (DL_FUNC) &_TreeDist_cpp_matching_split_distance, 3},
    {"_TreeDist_cpp_jaccard_similarity", (DL_FUNC) &_TreeDist_cpp_jaccard_similarity, 5},
    {"_TreeDist_cpp_msi_distance", (DL_FUNC) &_TreeDist_cpp_msi_distance, 3},
    {"_TreeDist_cpp_mutual_clustering", (DL_FUNC) &_TreeDist_cpp_mutual_clustering, 3},
    {"_TreeDist_cpp_shared_phylo", (DL_FUNC) &_TreeDist_cpp_shared_phylo, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_TreeDist(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
