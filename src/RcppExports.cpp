// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// cpp_robinson_foulds_distance
List cpp_robinson_foulds_distance(NumericMatrix x, NumericMatrix y, NumericVector nTip);
RcppExport SEXP _TreeDist_cpp_robinson_foulds_distance(SEXP xSEXP, SEXP ySEXP, SEXP nTipSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type nTip(nTipSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_robinson_foulds_distance(x, y, nTip));
    return rcpp_result_gen;
END_RCPP
}
// cpp_matching_split_distance
List cpp_matching_split_distance(NumericMatrix x, NumericMatrix y, NumericVector nTip);
RcppExport SEXP _TreeDist_cpp_matching_split_distance(SEXP xSEXP, SEXP ySEXP, SEXP nTipSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type nTip(nTipSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_matching_split_distance(x, y, nTip));
    return rcpp_result_gen;
END_RCPP
}
// cpp_jaccard_distance
List cpp_jaccard_distance(NumericMatrix x, NumericMatrix y, NumericVector nTip, NumericVector k, LogicalVector arboreal);
RcppExport SEXP _TreeDist_cpp_jaccard_distance(SEXP xSEXP, SEXP ySEXP, SEXP nTipSEXP, SEXP kSEXP, SEXP arborealSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type nTip(nTipSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type k(kSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type arboreal(arborealSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_jaccard_distance(x, y, nTip, k, arboreal));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_TreeDist_cpp_robinson_foulds_distance", (DL_FUNC) &_TreeDist_cpp_robinson_foulds_distance, 3},
    {"_TreeDist_cpp_matching_split_distance", (DL_FUNC) &_TreeDist_cpp_matching_split_distance, 3},
    {"_TreeDist_cpp_jaccard_distance", (DL_FUNC) &_TreeDist_cpp_jaccard_distance, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_TreeDist(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
