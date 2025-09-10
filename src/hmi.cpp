#include "hpart.h"
#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

namespace TreeDist {

static inline double x_log_x(size_t x) {
  return x ? x * std::log(x) : 0.0;
}

static inline size_t intersection_size(const std::vector<uint64_t>& A,
                                       const std::vector<uint64_t>& B) {
  size_t count = 0;
  ASSERT(A.size() == B.size());
  size_t size = A.size();
  
  for (size_t i = 0; i < size; ++i) {
    count += __builtin_popcountll(A[i] & B[i]);
  }
  return count;
}

// Hierarchical Mutual Information core
std::pair<size_t, double> hierarchical_mutual_info(
    const std::vector<TreeDist::HNode>& u_nodes,
    size_t u_idx,
    const std::vector<TreeDist::HNode>& v_nodes,
    size_t v_idx
) {
  const auto& Ut = u_nodes[u_idx];
  const auto& Us = v_nodes[v_idx];
  
  if (Ut.all_kids_leaves || Us.all_kids_leaves) {
    return {intersection_size(Ut.bitset, Us.bitset), 0.0};
  }
  
  size_t n_ts = 0;
  double H_uv = 0.0;
  double H_us = 0.0;
  double H_tv = 0.0;
  const size_t Us_size = Us.children.size();
  std::vector<size_t> n_tv(Us_size, 0);
  double mean_I_ts = 0.0;
  
  for (size_t u_child_idx : Ut.children) {
    size_t n_us = 0;
    for (size_t v = 0; v < Us_size; ++v) {
      size_t v_child_idx = Us.children[v];
      auto [n_uv, I_uv] = hierarchical_mutual_info(u_nodes, u_child_idx,
                                                   v_nodes, v_child_idx);
      
      n_ts += n_uv;
      n_tv[v] += n_uv;
      n_us += n_uv;
      H_uv += x_log_x(n_uv);
      mean_I_ts += n_uv * I_uv;
    }
    H_us += x_log_x(n_us);
  }
  
  for (auto _n_tv : n_tv) {
    H_tv += x_log_x(_n_tv);
  }
  
  if (n_ts == 0) {
    return {0, 0.0};
  }
  
  double local_I_ts = std::log(static_cast<double>(n_ts)) - (H_us + H_tv - H_uv) / static_cast<double>(n_ts);
  mean_I_ts /= static_cast<double>(n_ts);
  double I_ts = local_I_ts + mean_I_ts;
  
  return {n_ts, I_ts};
}
double hierarchical_self_info(const std::vector<TreeDist::HNode>& nodes, size_t idx) {
  const auto& n = nodes[idx];
  
  if (n.all_kids_leaves) return 0.0;
  
  size_t n_ts = 0;
  double H_uv = 0.0;
  double H_us = 0.0;
  double H_tv = 0.0;
  const size_t n_children = n.children.size();
  std::vector<size_t> n_tv(n_children, 0);
  double mean_I_ts = 0.0;
  
  for (size_t i = 0; i < n_children; ++i) {
    size_t u_idx = n.children[i];
    size_t n_us = 0;
    for (size_t j = 0; j < n_children; ++j) {
      size_t v_idx = n.children[j];
      size_t n_uv = (u_idx == v_idx) ? nodes[u_idx].leaf_count
      : intersection_size(nodes[u_idx].bitset, nodes[v_idx].bitset);
      
      n_ts += n_uv;
      n_tv[j] += n_uv;
      n_us += n_uv;
      H_uv += x_log_x(n_uv);
      mean_I_ts += n_uv * hierarchical_self_info(nodes, u_idx); // recurse
    }
    H_us += x_log_x(n_us);
  }
  
  for (auto _n_tv : n_tv) H_tv += x_log_x(_n_tv);
  
  if (n_ts == 0) return 0.0;
  double local_I_ts = std::log(static_cast<double>(n_ts)) - (H_us + H_tv - H_uv) / static_cast<double>(n_ts);
  mean_I_ts /= static_cast<double>(n_ts);
  
  return local_I_ts + mean_I_ts;
}



} // namespace TreeDist

// [[Rcpp::export]]
double HMI_xptr(SEXP ptr1, SEXP ptr2) {
  Rcpp::XPtr<TreeDist::HPart> hp1(ptr1);
  Rcpp::XPtr<TreeDist::HPart> hp2(ptr2);
  return TreeDist::hierarchical_mutual_info(hp1->nodes, hp1->root,
                                            hp2->nodes, hp2->root).second;
}

// [[Rcpp::export]]
double HME_xptr(SEXP ptr) {
  Rcpp::XPtr<TreeDist::HPart> hp(ptr);
  return hierarchical_self_info(hp->nodes, hp->root);
}
