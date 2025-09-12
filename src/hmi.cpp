#include "hpart.h"
#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

namespace TreeDist {

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
  
  const size_t Us_size = Us.children.size();
  
  size_t n_ts = 0;
  double H_uv = 0.0;
  double H_us = 0.0;
  double H_tv = 0.0;
  std::vector<size_t> n_tv(Us_size, 0);
  double mean_I_ts = 0.0;
  
  for (const size_t u_child_idx : Ut.children) {
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
  
  for (const auto _n_tv : n_tv) {
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

double character_mutual_info(
    const std::vector<TreeDist::HNode>& nodes, const size_t idx,
    const std::vector<std::vector<uint64_t>>& bitsets
) {
  
  const auto& nd = nodes[idx];
  if (nd.leaf_count < 2) return 0; // Exit early
  
  if (nd.all_kids_leaves) {
    const auto& nd_bits = nd.bitset;
    double h = nd.x_log_x;
    for (const auto& chr_bits : bitsets) {
      const size_t n = intersection_size(nd_bits, chr_bits);
      // Rcpp::Rcout << " Cherry " << idx << ": Intersection of " <<
      //   " (" << intersection_size(nd_bits, nd_bits) << ", " << 
      //     intersection_size(chr_bits, chr_bits) << ") = " << n << ".\n";
      
      // 1b. Continue sum of node's joint information
      h -= x_log_x(n);
    }
    // Rcpp::Rcout << "    Cherry " << idx << ", h = " << (h / std::log(2)) << ".\n";
    return h;
  }
  
  // Joint info = n_tips * sum [x_log_x(confusion_matrix_tips / n_tips)]
  //            = x_log_x(n_tips) - sum [x_log_x(confusion_matrix_tips)]
  double h = nd.x_log_x;

  // Rcpp::Rcout << "\n Node " << idx << ": initialize H to " << nd.leaf_count <<
  //   " log2(" << nd.leaf_count << ") = " << (h / std::log(2)) << ". \n";
  
  for (const auto& child : nd.children) {
    const auto& cld_bits = nodes[child].bitset;
    // Rcpp::Rcout << " < Before child " << child << ", h = " << (h / std::log(2)) << ".\n";
    for (const auto& chr_bits : bitsets) {
      // 1a. Populate cell in confusion matrix
      const size_t n = intersection_size(cld_bits, chr_bits);
      // Rcpp::Rcout << "     Child " << child << ": Intersection of " <<
      //   " (" << intersection_size(cld_bits, cld_bits) << ", " << 
      //   intersection_size(chr_bits, chr_bits) << ") = " << n << ".\n";
          
      // 1b. Continue sum of node's joint information
      h -= x_log_x(n);
    }
    // Rcpp::Rcout << " > After child " << child << ", h = " << (h / std::log(2)) << ".\n";
  }
  
  // Rcpp::Rcout << " Unconditioned joint H(char, " << idx << ") = " <<
  //   (h / std::log(2)) << " bits.\n\n";
  
  
  
  for (const auto& child : nd.children) { // TODO reintegrate into single loop
    // 2. Load the contributions to tree entropy from child nodes, in postorder.
    const auto& cld_bits = nodes[child].bitset;
    
    if (nodes[child].leaf_count > 1) { // TODO: can we revert to child_h > 0?
      
      // Remove joint info we've already counted in the parent:
      h -= nodes[child].x_log_x;
      for (const auto& chr_bits : bitsets) {
        const size_t n = intersection_size(cld_bits, chr_bits);
        // Rcpp::Rcout << "     Child " << child << ": Intersection of " <<
        //   " (" << intersection_size(cld_bits, cld_bits) << ", " << 
        //     intersection_size(chr_bits, chr_bits) << ") = " << n << 
        //       "; reducing H by " << x_log_x(n) << ".\n";
        h += x_log_x(n);
      }
      
      // Then add the unconditioned information contained within the child
      // subtree
      // 
      // Propagate in postorder
      const double child_contribuition = character_mutual_info(nodes, child, bitsets);
      // Rcpp::Rcout << " Adding subtree contribution from " << child << " = " <<
      //   (child_contribuition / std::log(2))<< "\n";
      h += child_contribuition;
    }
  }
  
  // Rcpp::Rcout << " >>> Final h below " << idx << " is " << (h / std::log(2)) << ".\n\n";
  return h;
}

double tree_entropy(const std::vector<TreeDist::HNode>& nodes,
                    const size_t idx) {
  const auto& nd = nodes[idx];
  
  if (nd.all_kids_leaves) return 0; // As we don't distinguish these leaves
  // We might expect to return nd.x_log_x: but consider the star tree,
  // whose information content we expect to be zero.
  // We don't think of the star tree as assigning its leaves to labelled
  // size-1 clusters.
  
  // Info = n_tips * sum [x_log_x(child_i_tips / n_tips)]
  //      = x_log_x(n_tips) - sum [x_log_x(child_i_tips)] <- easier to calculate
  
  double h = nd.x_log_x;
  
  for (const auto& child : nd.children) {
    
    // 1. Continue the summation of this node's entropy
    const double child_h = nodes[child].x_log_x;
    h -= child_h; // contribution to this node's entropy
    
    // 2. Load the contributions to tree entropy from child nodes, in postorder.
    if (child_h > 0) {
      const double child_contribuition = tree_entropy(nodes, child);
      h += child_contribuition;
    }
  }
  
  return h;
}



} // namespace TreeDist

// [[Rcpp::export]]
double HMI_xptr(SEXP ptr1, SEXP ptr2) {
  Rcpp::XPtr<TreeDist::HPart> hp1(ptr1);
  Rcpp::XPtr<TreeDist::HPart> hp2(ptr2);
  if (hp1->nodes[hp1->root].n_tip != hp2->nodes[hp2->root].n_tip) {
    Rcpp::stop("Trees must have the same number of leaves");
  }
  return TreeDist::hierarchical_mutual_info(hp1->nodes, hp1->root,
                                            hp2->nodes, hp2->root).second;
}

// [[Rcpp::export]]
double HH_xptr(SEXP ptr) {
  Rcpp::XPtr<TreeDist::HPart> hp(ptr);
  if (hp->hier_entropy == std::numeric_limits<double>::min()) {
    // When requring C++26, update to constexpr
    const double eps = std::sqrt(std::numeric_limits<double>::epsilon());
    const double value = hierarchical_self_info(hp->nodes, hp->root);
    hp->hier_entropy = std::abs(value) < eps ? 0 : value;
  }
  return hp->hier_entropy;
}

//' Directly calculate entropy from an HPart pointer
//' 
//' Intended for developer use only; no safeguards against bad input.
//' May crash R if used incorrectly.
//' 
//' @param ptr,char_ptr,tree_ptr Pointers to HPart objects for which chosen
//' entropy value should be calculated.
//' @inheritParams CharAMI
//' 
//' @template MRS
//' @keywords internal
//' @export
// [[Rcpp::export]]
double H_xptr(SEXP ptr) {
  Rcpp::XPtr<TreeDist::HPart> hp(ptr);
  if (hp->entropy == std::numeric_limits<double>::min()) {
    // When requring C++26, update to constexpr
    const double eps = std::sqrt(std::numeric_limits<double>::epsilon());
    const double value = tree_entropy(hp->nodes, hp->root);
    hp->entropy = std::abs(value) < eps ? 0 : value;
  }
  return hp->entropy;
}

//' @rdname H_xptr
//' @export
// [[Rcpp::export]]
double JH_xptr(SEXP char_ptr, SEXP tree_ptr) {
  Rcpp::XPtr<TreeDist::HPart> ch(char_ptr);
  Rcpp::XPtr<TreeDist::HPart> tr(tree_ptr);
  if (ch->nodes[ch->root].n_tip != tr->nodes[tr->root].n_tip) {
    Rcpp::stop("Tree and character must describe the same leaves");
  }
  std::vector<std::vector<uint64_t>> bitsets;
  auto ch_root = ch->nodes[ch->root];
  bitsets.reserve(ch_root.children.size());
  bool warned = false;
  for (const auto& state : ch_root.children) {
    bitsets.push_back(ch->nodes[state].bitset);
    if (!warned && !ch->nodes[state].all_kids_leaves) {
      Rcpp::warning("Character is a tree; only first level of hierarchy used");
      warned = true;
    }
  }
  
  return TreeDist::character_mutual_info(tr->nodes, tr->root, bitsets);
}

inline void fisher_yates_shuffle(std::vector<int>& v) noexcept {
  for (size_t i = v.size() - 1; i > 0; --i) {
    size_t j = static_cast<size_t>(std::floor(R::unif_rand() * (i + 1)));
    std::swap(v[i], v[j]);
  }
}

// [[Rcpp::export]]
Rcpp::NumericVector EHMI_xptr(SEXP hp1_ptr, SEXP hp2_ptr,
                              double tolerance = 0.01,
                              int minResample = 36) {
  
  if (minResample < 2) {
    Rcpp::stop("Must perform at least one resampling");
  }
  if (tolerance < 1e-8) {
    Rcpp::stop("Tolerance too low");
  }
  
  Rcpp::XPtr<TreeDist::HPart> hp1(hp1_ptr);
  Rcpp::XPtr<TreeDist::HPart> hp2(hp2_ptr);
  
  const size_t n_tip = hp1->nodes[hp1->root].n_tip;
  ASSERT(hp2->nodes[hp2->root].n_tip == n_tip);
  
  double runMean = 0.0;
  double runS = 0.0;
  int runN = 0;
  double relativeError = tolerance * 2; // Avoid -Wmaybe-uninitialized
  
  Rcpp::RNGScope scope;
  
  SEXP hp1_shuf = clone_hpart(hp1_ptr);
  std::vector<int> shuffled(n_tip);
  std::iota(shuffled.begin(), shuffled.end(), 0);
  
  while (relativeError > tolerance || runN < minResample) {
    // Shuffle leaves
    fisher_yates_shuffle(shuffled);
    
    // Apply shuffled labels
    relabel_hpart(hp1_shuf, shuffled);
    
    // Compute HMI
    double x = HMI_xptr(hp1_shuf, hp2);
    
    // Welford update
    runN++;
    double delta = x - runMean;
    runMean += delta / runN;
    runS += delta * (x - runMean);
    
    double runVar = (runN > 1) ? runS / (runN - 1) : 0.0;
    double runSD = std::sqrt(runVar);
    double runSEM = runSD / std::sqrt(runN);
    relativeError = std::abs(runMean) < 1e-6 ?
      runSEM :
      runSEM / std::abs(runMean);
  }
  
  double runVar = (runN > 1) ? runS / (runN - 1) : 0.0;
  double runSD = std::sqrt(runVar);
  double runSEM = runSD / std::sqrt(runN);
  
  Rcpp::NumericVector result = Rcpp::NumericVector::create(runMean);
  result.attr("var") = runVar;
  result.attr("sd") = runSD;
  result.attr("sem") = runSEM;
  result.attr("samples") = runN;
  result.attr("relativeError") = relativeError;
  
  return result;
}

//' @rdname H_xptr
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector EJH_xptr(SEXP char_ptr, SEXP tree_ptr,
                             double tolerance = 0.01,
                             int minResample = 36) {
  
  if (minResample < 2) {
    Rcpp::stop("Must perform at least one resampling");
  }
  if (tolerance < 1e-8) {
    Rcpp::stop("Tolerance too low");
  }
  
  Rcpp::XPtr<TreeDist::HPart> ch(char_ptr);
  Rcpp::XPtr<TreeDist::HPart> tr(tree_ptr);
  
  const size_t n_tip = ch->nodes[ch->root].n_tip;
  if (tr->nodes[tr->root].n_tip != static_cast<int>(n_tip)) {
    Rcpp::stop("Tree and character must describe the same leaves");
  }
  
  double runMean = 0.0;
  double runS = 0.0;
  int runN = 0;
  double relativeError = tolerance * 2; // Avoid -Wmaybe-uninitialized
  
  Rcpp::RNGScope scope;
  
  SEXP ch_shuf = clone_hpart(char_ptr);
  std::vector<int> shuffled(n_tip);
  std::iota(shuffled.begin(), shuffled.end(), 0);
  
  while (relativeError > tolerance || runN < minResample) {
    // Shuffle leaves
    fisher_yates_shuffle(shuffled);
    
    // Apply shuffled labels
    relabel_hpart(ch_shuf, shuffled);
    
    // Compute JH
    double x = JH_xptr(ch_shuf, tr);
    
    // Welford update
    runN++;
    double delta = x - runMean;
    runMean += delta / runN;
    runS += delta * (x - runMean);
    
    double runVar = (runN > 1) ? runS / (runN - 1) : 0.0;
    double runSD = std::sqrt(runVar);
    double runSEM = runSD / std::sqrt(runN);
    relativeError = std::abs(runMean) < 1e-6 ?
      runSEM :
      runSEM / std::abs(runMean);
  }
  
  double runVar = (runN > 1) ? runS / (runN - 1) : 0.0;
  double runSD = std::sqrt(runVar);
  double runSEM = runSD / std::sqrt(runN);
  
  Rcpp::NumericVector result = Rcpp::NumericVector::create(runMean);
  result.attr("var") = runVar;
  result.attr("sd") = runSD;
  result.attr("sem") = runSEM;
  result.attr("samples") = runN;
  result.attr("relativeError") = relativeError;
  
  return result;
}
