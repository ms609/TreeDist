// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// COMCLUST
IntegerMatrix COMCLUST(List trees);
RcppExport SEXP _TreeDist_COMCLUST(SEXP treesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type trees(treesSEXP);
    rcpp_result_gen = Rcpp::wrap(COMCLUST(trees));
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
// cpp_mmsi_distance
List cpp_mmsi_distance(const RawMatrix x, const RawMatrix y, const IntegerVector nTip);
RcppExport SEXP _TreeDist_cpp_mmsi_distance(SEXP xSEXP, SEXP ySEXP, SEXP nTipSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const RawMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< const RawMatrix >::type y(ySEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type nTip(nTipSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_mmsi_distance(x, y, nTip));
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
    {"_TreeDist_COMCLUST", (DL_FUNC) &_TreeDist_COMCLUST, 1},
    {"_TreeDist_lapjv", (DL_FUNC) &_TreeDist_lapjv, 2},
    {"_TreeDist_cpp_mast", (DL_FUNC) &_TreeDist_cpp_mast, 3},
    {"_TreeDist_cpp_nni_distance", (DL_FUNC) &_TreeDist_cpp_nni_distance, 3},
    {"_TreeDist_cpp_robinson_foulds_distance", (DL_FUNC) &_TreeDist_cpp_robinson_foulds_distance, 3},
    {"_TreeDist_cpp_robinson_foulds_info", (DL_FUNC) &_TreeDist_cpp_robinson_foulds_info, 3},
    {"_TreeDist_cpp_matching_split_distance", (DL_FUNC) &_TreeDist_cpp_matching_split_distance, 3},
    {"_TreeDist_cpp_jaccard_similarity", (DL_FUNC) &_TreeDist_cpp_jaccard_similarity, 5},
    {"_TreeDist_cpp_mmsi_distance", (DL_FUNC) &_TreeDist_cpp_mmsi_distance, 3},
    {"_TreeDist_cpp_mutual_clustering", (DL_FUNC) &_TreeDist_cpp_mutual_clustering, 3},
    {"_TreeDist_cpp_shared_phylo", (DL_FUNC) &_TreeDist_cpp_shared_phylo, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_TreeDist(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
