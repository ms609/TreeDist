#include "hpart.h"
#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// Compute entropy of a block at node i in hp
static inline double entropy_node(const Node &nd, size_t nTips) {
  double p = static_cast<double>(nd.leafCount) / nTips;
  return -p * std::log(p);
}
// A helper function to compute popcount of the intersection
static inline size_t intersection_size(
    const std::vector<uint64_t>& bitset1,
    const std::vector<uint64_t>& bitset2) {
  size_t count = 0;
  ASSERT(bitset1.size() == bitset2.size());
  size_t size = bitset1.size();
  
  for (size_t i = 0; i < size; ++i) {
    uint64_t intersection_bits = bitset1[i] & bitset2[i];
    count += __builtin_popcountll(intersection_bits);
  }
  return count;
}

const inline double x_log_x(size_t x) {
  return x ? x * std::log(x) : 0.0;
}

// Hierarchical Mutual Information core
std::pair<size_t, double> hierarchical_mutual_info(const Node &Ut, const Node &Us) {
  if (Ut.allKidsLeaves || Us.allKidsLeaves) {
    return std::pair<size_t, double>{intersection_size(Ut.bitset, Us.bitset),
                                     0.0};
  }
  size_t n_ts = 0;
  double H_uv = 0.0;
  double H_us = 0.0;
  double H_tv = 0.0;
  const size_t Us_size = Us.children.size();
  std::vector<size_t> n_tv(Us_size);
  double mean_I_ts = 0.0;
  for (const auto& Uu : Ut.children) {
    size_t n_us = 0.0;
    for (size_t v = 0; v < Us_size; ++v) {
      const auto& Uv = Us.children[v];
      const std::pair<size_t, double> niUV = hierarchical_mutual_info(Uu, Uv);
      const size_t n_uv = niUV.first;
      const double I_uv = niUV.second;
      n_ts += n_uv;
      n_tv[v] += n_uv;
      n_us += n_uv;
      H_uv += x_log_x(n_uv);
      mean_I_ts += n_uv * I_uv;
    }
    H_us += x_log_x(n_us);
  }
  for (const auto& _n_tv : n_tv) {
    H_tv += x_log_x(_n_tv);
  }
  if (n_ts == 0) {
    return std::pair<size_t, double>{0, 0.0};
  }
  const double local_I_ts = std::log(n_ts) - (H_us + H_tv - H_uv) /
    static_cast<double>(n_ts);
  mean_I_ts /= static_cast<double>(n_ts);
  const double I_ts = local_I_ts + mean_I_ts;
  return std::pair<size_t, double>{n_ts, I_ts};
}

// [[Rcpp::export]]
double HMI_xptr(SEXP ptr1, SEXP ptr2) {
  Rcpp::XPtr<HPart> hp1(ptr1);
  Rcpp::XPtr<HPart> hp2(ptr2);
  return hierarchical_mutual_info(hp1.nodes[hp1.rootIndex],
                                  hp2.nodes[hp2.rootIndex]).second;
}
