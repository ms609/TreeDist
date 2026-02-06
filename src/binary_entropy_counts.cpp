// src/binary_entropy_counts.cpp
#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// Compute H(p) = -p log p - (1-p) log(1-p) for p = count / nLeaves
// Counts at 0 or nLeaves get entropy 0; NA counts return NA.
// [[Rcpp::export]]
NumericVector binary_entropy_counts(IntegerVector inSplit, int nLeaves) {
  int K = inSplit.size();
  NumericVector out(K);
  
  if (nLeaves <= 0) {
    std::fill(out.begin(), out.end(), NA_REAL);
    return out;
  }
  
  const double invN = 1.0 / static_cast<double>(nLeaves);
  
  for (int j = 0; j < K; ++j) {
    int a = inSplit[j];
    if (IntegerVector::is_na(a)) {
      out[j] = NA_REAL;
      continue;
    }
    if (a <= 0 || a >= nLeaves) {
      // p ∈ {0,1} ⇒ H = 0
      out[j] = 0.0;
      continue;
    }
    
    double x = a * invN;
    const double one_minus_x = 1.0 - x;
    
    // Use log1p(-x) for better precision near x≈1
    double H = -(x * std::log(x) + one_minus_x * std::log1p(-x));
    out[j] = H;
  }
  constexpr double oneOverLog2 = 1.442695040888963387005;
  return out * oneOverLog2;
}
