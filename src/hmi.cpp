#include "hpart.h"
#include <Rcpp.h>
#include <cmath>
#include <algorithm> // for std::shuffle
#include <RcppParallel.h>
#include <random>
#include <thread>
using namespace Rcpp;

namespace TreeDist {

static inline size_t intersection_size(const uint64_t* A,
                                       const uint64_t* B,
                                       const size_t size) {
  return std::transform_reduce(
    A, A + size, B,
    size_t{0},      // initial value
    std::plus<>{},  // reduction operation
    [](uint64_t a, uint64_t b) { return __builtin_popcountll(a & b); }  // transform
  );
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
    return {intersection_size(Ut.bitset, Us.bitset, Us.n_block), 0.0};
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
      : intersection_size(nodes[u_idx].bitset, nodes[v_idx].bitset,
        nodes[u_idx].n_block);
      
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
    const std::vector<const uint64_t*>& bitsets
) {
  
  const auto& nd = nodes[idx];
  if (nd.leaf_count < 2) return 0; // Exit early
  const size_t n_block = nodes[idx].n_block;
  
  if (nd.all_kids_leaves) {
    const uint64_t* nd_bits = nd.bitset;
    double h = nd.x_log_x;
    for (const uint64_t* chr_bits : bitsets) {
      const size_t n = intersection_size(nd_bits, chr_bits, n_block);
      h -= x_log_x(n);
    }
    return h;
  }
  
  // Joint info = n_tips * sum [x_log_x(confusion_matrix_tips / n_tips)]
  //            = x_log_x(n_tips) - sum [x_log_x(confusion_matrix_tips)]
  double h = nd.x_log_x;
  
  for (const auto& child : nd.children) {
    const double child_h = nodes[child].x_log_x;
    if (child_h > 0) {
      h -= child_h;
      
      // We can simply add the information contained within the child subtree
      // The full computation of this node's entropy cancels with the
      // conditional entropy of the descendent node.
      // 
      // Propagate in postorder
      const double child_contribution = character_mutual_info(nodes, child,
                                                              bitsets);
      h += child_contribution;
    }
  }
  
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

inline double JH_core(const TreeDist::HPart* ch,
                      const TreeDist::HPart* tr) {
  
  const auto& ch_root = ch->nodes[ch->root];
  std::vector<const uint64_t*> bitsets;
  bitsets.reserve(ch_root.children.size());
  
  for (const auto& state : ch_root.children) {
    bitsets.push_back(ch->nodes[state].bitset);
  }
  
  return TreeDist::character_mutual_info(tr->nodes, tr->root, bitsets);
}


//' @rdname H_xptr
//' @export
// [[Rcpp::export]]
double JH_xptr(SEXP char_ptr, SEXP tree_ptr) {
  Rcpp::XPtr<TreeDist::HPart> ch(char_ptr);
  Rcpp::XPtr<TreeDist::HPart> tr(tree_ptr);
  if (ch->nodes[ch->root].n_tip != tr->nodes[tr->root].n_tip) {
    Rcpp::Rcout << "Character: " << ch->nodes[ch->root].n_tip << 
      "; tree: " << tr->nodes[tr->root].n_tip << "|\n";
    Rcpp::stop("Tree and character must describe the same leaves");
  }
  auto ch_root = ch->nodes[ch->root];
  bool warned = false;
  for (const auto& state : ch_root.children) {
    if (!warned && !ch->nodes[state].all_kids_leaves) {
      Rcpp::warning("Character is a tree; only first level of hierarchy used");
      warned = true;
    }
  }
  
  return JH_core(ch.get(), tr.get());
}

// [[Rcpp::export]]
Rcpp::NumericVector EHMI_xptr(const SEXP hp1_ptr, const SEXP hp2_ptr,
                              const double precision = 0.01,
                              const int minResample = 36) {
  
  if (minResample < 2) {
    Rcpp::stop("Must perform at least one resampling");
  }
  if (precision < 1e-8) {
    Rcpp::stop("Precision too small");
  }
  
  Rcpp::XPtr<TreeDist::HPart> hp1(hp1_ptr);
  Rcpp::XPtr<TreeDist::HPart> hp2(hp2_ptr);
  
  const size_t n_tip = hp1->nodes[hp1->root].n_tip;
  ASSERT(hp2->nodes[hp2->root].n_tip == n_tip);
  
  double run_mean = 0.0;
  double run_s = 0.0;
  int run_n = 0;
  double relative_error = precision * 2; // Avoid -Wmaybe-uninitialized
  
  SEXP hp1_shuf = clone_hpart(hp1_ptr);
  std::vector<int> shuffled(n_tip);
  std::iota(shuffled.begin(), shuffled.end(), 0);
  
  Rcpp::RNGScope scope;
  unsigned int seed =
    static_cast<unsigned int>(R::unif_rand() * 
    std::numeric_limits<unsigned int>::max());
  std::mt19937_64 rng(seed);
  
  while (relative_error > precision || run_n < minResample) {
    
    std::shuffle(shuffled.begin(), shuffled.end(), rng);
    relabel_hpart(hp1_shuf, shuffled);
    
    double x = HMI_xptr(hp1_shuf, hp2);
    
    // Welford update
    run_n++;
    double delta = x - run_mean;
    run_mean += delta / run_n;
    run_s += delta * (x - run_mean);
    
    double run_var = (run_n > 1) ? run_s / (run_n - 1) : 0.0;
    double run_sd = std::sqrt(run_var);
    double run_sem = run_sd / std::sqrt(run_n);
    relative_error = std::abs(run_mean) < 1e-6 ?
      run_sem :
      run_sem / std::abs(run_mean);
  }
  
  double run_var = (run_n > 1) ? run_s / (run_n - 1) : 0.0;
  double run_sd = std::sqrt(run_var);
  double run_sem = run_sd / std::sqrt(run_n);
  
  Rcpp::NumericVector result = Rcpp::NumericVector::create(run_mean);
  result.attr("var") = run_var;
  result.attr("sd") = run_sd;
  result.attr("sem") = run_sem;
  result.attr("samples") = run_n;
  result.attr("precision") = relative_error;
  
  return result;
}

#include <RcppParallel.h>
#include <random>

struct MCChain : public RcppParallel::Worker {
  
  // Input data
  const TreeDist::HPart* ch;
  const TreeDist::HPart* tr;
  const double precision;
  const int minResample;
  const std::function<double(const double,const double)> propagate_sem;
  const std::function<double(const double,const double)> rel_err;
  const unsigned int base_seed;
  const size_t n_tip;
  
  // Output data (one per thread)
  RcppParallel::RVector<double> results_mean;
  RcppParallel::RVector<double> results_var;
  RcppParallel::RVector<int> results_n;
  
  MCChain(const TreeDist::HPart* ch, const TreeDist::HPart* tr,
          const double precision,
          const int minResample,
          const std::function<double(const double,const double)>& propagate_sem,
          const std::function<double(const double,const double)>& rel_err,
          const unsigned int base_seed, size_t n_tip,
          Rcpp::NumericVector results_mean,
          Rcpp::NumericVector results_var,
          Rcpp::IntegerVector results_n) :
    ch(ch), tr(tr), precision(precision), minResample(minResample),
    propagate_sem(propagate_sem), rel_err(rel_err),
    base_seed(base_seed), n_tip(n_tip), results_mean(results_mean),
    results_var(results_var), results_n(results_n) {}
  
  void operator()(std::size_t begin, std::size_t end) {
    
    // Each worker processes one chain (begin to end should be consecutive indices)
    for (std::size_t chain_id = begin; chain_id < end; chain_id++) {
      
      // Thread-specific RNG
      std::mt19937_64 rng(base_seed + static_cast<unsigned int>(chain_id));
      
      // Thread-specific character copy and shuffle vector
      
      TreeDist::HPart* ch_shuf = new TreeDist::HPart(*ch);
      std::vector<int> shuffled(n_tip);
      std::iota(shuffled.begin(), shuffled.end(), 0);
      
      // Chain statistics
      double run_mean = 0.0;
      double run_s = 0.0;
      int run_n = 0;
      double relative_error = precision * 2;
      
      double run_var, run_sd, jh_sem, run_sem;
      
      // Monte Carlo loop (same logic as original)
      while (relative_error > precision || run_n < minResample) {
        
        std::shuffle(shuffled.begin(), shuffled.end(), rng);
        TreeDist::relabel_hpart_xptr(ch_shuf, shuffled);
        
        double x = JH_core(ch_shuf, tr);
        
        // Welford update
        run_n++;
        double delta = x - run_mean;
        run_mean += delta / run_n;
        run_s += delta * (x - run_mean);
        
        run_var = (run_n > 1) ? run_s / (run_n - 1) : 0.0;
        run_sd = std::sqrt(run_var);
        jh_sem = run_sd / std::sqrt(run_n);
        run_sem = propagate_sem(run_mean, jh_sem);
        relative_error = rel_err(run_mean, run_sem);
      }
      
      // Store results
      results_mean[chain_id] = run_mean;
      results_var[chain_id] = run_var;
      results_n[chain_id] = run_n;
    }
  }
};

Rcpp::NumericVector EJH_core_parallel(const TreeDist::HPart* ch,
                                      const TreeDist::HPart* tr,
                                      const unsigned int base_seed,
                                      const double precision = 0.01,
                                      const int minResample = 36,
                                      int n_threads = -1,
                                      std::function<double(const double,const double)> propagate_sem =
                                        [](const double jh_mean, const double jh_sem) {
                                          return jh_sem;
                                        },
                                        std::function<double(const double,const double)> rel_err =
                                          [](const double jh_mean, const double run_sem) {
                                            return std::abs(jh_mean) < 1e-6 ?
                                            run_sem :
                                            run_sem / std::abs(jh_mean);
                                          }
) {
  
  const size_t n_tip = ch->nodes[ch->root].n_tip;
  if (tr->nodes[tr->root].n_tip != static_cast<int>(n_tip)) {
    Rcpp::stop("Tree and character must describe the same leaves");
  }
  
  // Determine number of threads
  if (n_threads <= 0) {
    n_threads = std::thread::hardware_concurrency();
  }
  if (n_threads > (int)std::thread::hardware_concurrency()) {
    n_threads = std::thread::hardware_concurrency();
  }
  
  // Scale precision for parallel chains
  double chain_precision = precision * std::sqrt(static_cast<double>(n_threads));
  
  // Storage for chain results
  Rcpp::NumericVector chain_means(n_threads);
  Rcpp::NumericVector chain_vars(n_threads);
  Rcpp::IntegerVector chain_ns(n_threads);
  
  // Create and run parallel worker
  MCChain worker(ch, tr, chain_precision, minResample, propagate_sem, rel_err,
                 base_seed, n_tip, chain_means, chain_vars, chain_ns);
  
  RcppParallel::parallelFor(0, n_threads, worker);
  
  // Combine results from all chains
  double combined_mean = 0.0;
  double combined_var = 0.0;
  int total_samples = 0;
  
  // Simple average of means
  for (int i = 0; i < n_threads; i++) {
    combined_mean += chain_means[i];
    total_samples += chain_ns[i];
  }
  combined_mean /= n_threads;
  
  // Combined variance: Var[X̄] = (1/n²) * Σ Var[X_i] for independent X_i
  for (int i = 0; i < n_threads; i++) {
    combined_var += chain_vars[i] / (chain_ns[i] * n_threads * n_threads);
  }
  
  double combined_sd = std::sqrt(combined_var * total_samples);
  double combined_sem = combined_sd / std::sqrt(total_samples);
  double final_sem = propagate_sem(combined_mean, combined_sem);
  double final_precision = rel_err(combined_mean, final_sem);
  
  // Return result in same format as original
  Rcpp::NumericVector result = Rcpp::NumericVector::create(combined_mean);
  result.attr("samples") = total_samples;
  result.attr("ejh") = combined_mean;
  result.attr("ejhVar") = combined_var * total_samples;  // Scale back to total variance
  result.attr("ejhSD") = combined_sd;
  result.attr("ejhSEM") = combined_sem;
  result.attr("sem") = final_sem;
  result.attr("precision") = final_precision;
  result.attr("chains") = n_threads;
  
  return result;
}

Rcpp::NumericVector EJH_core(const TreeDist::HPart* ch,
                             const TreeDist::HPart* tr,
                             const unsigned int seed,
                             const double precision = 0.01,
                             const int minResample = 36,
                             std::function<double(const double,const double)> propagate_sem =
                               [](const double jh_mean, const double jh_sem) {
                                 return jh_sem;
                               },
                              std::function<double(const double,const double)> rel_err =
                               [](const double jh_mean, const double run_sem) {
                                 return std::abs(jh_mean) < 1e-6 ?
                                 run_sem :
                                 run_sem / std::abs(jh_mean);
                               }
                               ) {
  
  const size_t n_tip = ch->nodes[ch->root].n_tip;
  
  double run_mean = 0.0;
  double run_s = 0.0;
  int run_n = 0;
  double relative_error = precision * 2; // Avoid -Wmaybe-uninitialized
  
  TreeDist::HPart* ch_shuf = new TreeDist::HPart(*ch);
  std::vector<int> shuffled(n_tip);
  std::iota(shuffled.begin(), shuffled.end(), 0);
  
  std::mt19937_64 rng(seed);
  
  double run_var;
  double run_sd;
  double jh_sem;
  double run_sem;
  
  while (relative_error > precision || run_n < minResample) {
    
    std::shuffle(shuffled.begin(), shuffled.end(), rng);
    TreeDist::relabel_hpart_xptr(ch_shuf, shuffled);
    
    double x = JH_core(ch_shuf, tr);
    
    // Welford update
    run_n++;
    double delta = x - run_mean;
    run_mean += delta / run_n;
    run_s += delta * (x - run_mean);
    
    run_var = (run_n > 1) ? run_s / (run_n - 1) : 0.0;
    run_sd = std::sqrt(run_var);
    jh_sem = run_sd / std::sqrt(run_n);
    run_sem = propagate_sem(run_mean, jh_sem);
    relative_error = rel_err(run_mean, run_sem);
  }
  
  
  Rcpp::NumericVector result = Rcpp::NumericVector::create(run_mean);
  result.attr("samples") = run_n;
  result.attr("ejh") = run_mean;
  result.attr("ejhVar") = run_var;
  result.attr("ejhSD") = run_sd;
  result.attr("ejhSEM") = jh_sem;
  result.attr("sem") = run_sem;
  result.attr("precision") = relative_error;
  
  return result;
}

//' @rdname H_xptr
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector EJH_xptr(const SEXP char_ptr, const SEXP tree_ptr,
                             const double precision = 0.01,
                             const int minResample = 36,
                             const int nCores = 1) {
  // Input validation
  if (minResample < 2) {
    Rcpp::stop("Must perform at least one resampling");
  }
  if (precision < 1e-8) {
    Rcpp::stop("Precision too small");
  }
  
  Rcpp::XPtr<TreeDist::HPart> ch(char_ptr);
  Rcpp::XPtr<TreeDist::HPart> tr(tree_ptr);
  
  const size_t n_tip = ch->nodes[ch->root].n_tip;
  if (tr->nodes[tr->root].n_tip != static_cast<int>(n_tip)) {
    Rcpp::stop("Tree and character must describe the same leaves");
  }
  
  const unsigned int seed =
    static_cast<unsigned int>(R::unif_rand() * 
    std::numeric_limits<unsigned int>::max());
  
  return nCores < 2 ?
  EJH_core(ch.get(), tr.get(), seed, precision, minResample) :
  EJH_core_parallel(ch.get(), tr.get(), seed, precision, minResample, nCores);
}

//' @rdname H_xptr
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector EMI_xptr(const SEXP char_ptr, const SEXP tree_ptr,
                             const double precision = 0.01,
                             const int minResample = 36,
                             const int nCores = 1) {
  
  if (minResample < 2) {
    Rcpp::stop("Must perform at least one resampling");
  }
  if (precision < 1e-8) {
    Rcpp::stop("Precision too small");
  }
  
  Rcpp::XPtr<TreeDist::HPart> ch(char_ptr);
  Rcpp::XPtr<TreeDist::HPart> tr(tree_ptr);
  
  const size_t n_tip = ch->nodes[ch->root].n_tip;
  if (tr->nodes[tr->root].n_tip != static_cast<int>(n_tip)) {
    Rcpp::stop("Tree and character must describe the same leaves");
  }
  
  const double h1 = H_xptr(ch);
  const double h2 = H_xptr(tr);
  const double h1_h2 = h1 + h2;
  
  Rcpp::RNGScope scope;
  const unsigned int seed =
    static_cast<unsigned int>(R::unif_rand() * 
    std::numeric_limits<unsigned int>::max());
  
  auto emi_sem = [=](const double eh12, const double eh12_sem) {
    return eh12_sem;
  };
  
  auto emi_rel_err = [=](const double eh12, const double run_sem) {
    const double emi = h1_h2 - eh12;
    return std::abs(emi) < 1e-6 ?
      run_sem :
      run_sem / std::abs(emi);
  };
  
  
  Rcpp::NumericVector res = nCores < 2 ?
    EJH_core(ch.get(), tr.get(), seed, precision, minResample, emi_sem,
           emi_rel_err) :
    EJH_core_parallel(ch.get(), tr.get(), seed, precision, minResample, nCores,
                      emi_sem, emi_rel_err);
  
  const double emi = h1_h2 - res[0];
  res[0] = emi;
  
  return res;
}

//' @rdname H_xptr
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector AMI_xptr(const SEXP char_ptr, const SEXP tree_ptr,
                             const SEXP mean_fn, const double precision = 0.01,
                             int minResample = 36, const int nCores = 1) {
  
  if (minResample < 2) {
    Rcpp::stop("Must perform at least one resampling");
  }
  if (precision < 1e-8) {
    Rcpp::stop("Precision too small");
  }
  
  Rcpp::XPtr<TreeDist::HPart> ch(char_ptr);
  Rcpp::XPtr<TreeDist::HPart> tr(tree_ptr);
  
  const size_t n_tip = ch->nodes[ch->root].n_tip;
  if (tr->nodes[tr->root].n_tip != static_cast<int>(n_tip)) {
    Rcpp::stop("Tree and character must describe the same leaves");
  }
  
  const double h1 = H_xptr(ch);
  const double h2 = H_xptr(tr);
  
  Rcpp::Function mean(mean_fn);
  const double mn = Rcpp::as<double>(mean(h1, h2));
  const double h1_h2 = h1 + h2;
  
  const double h12 = JH_core(ch, tr);
  const double mi = h1_h2 - h12;
  
  const double eps = std::sqrt(std::numeric_limits<double>::epsilon());
  
  if (std::abs(mi - mn) < eps) {
    Rcpp::NumericVector result = Rcpp::NumericVector::create(1.0);
    result.attr("samples") = 0;
    result.attr("ejh") = NA_REAL;
    result.attr("ejhVar") = NA_REAL;
    result.attr("ejhSD") = NA_REAL;
    result.attr("ejhSEM") = NA_REAL;
    result.attr("sem") = 0;
    result.attr("precision") = 0;
    
    return result;
  }
  
  Rcpp::RNGScope scope;
  const unsigned int seed =
    static_cast<unsigned int>(R::unif_rand() * 
    std::numeric_limits<unsigned int>::max());

  auto ami_sem = [=](const double eh12, const double eh12_sem) {
    if (eh12_sem > eps) {
      const double emi = h1_h2 - eh12;
      const double deriv = (mi - mn) / ((mn - emi) * (mn - emi));
      const double ret = std::abs(deriv) * eh12_sem;
      return (ret < eps) ? 0.0 : ret;
    } else {
      return 0.0;
    }
  };
  
  auto absolute_err = [=](const double eh12, const double run_sem) {
    return run_sem;
  };
  
  
  Rcpp::NumericVector res = nCores < 2 ?
    EJH_core(ch.get(), tr.get(), seed, precision, minResample, ami_sem,
           absolute_err) :
    EJH_core_parallel(ch.get(), tr.get(), seed, precision, minResample, nCores,
                      ami_sem, absolute_err);
  
  const double emi = h1_h2 - res[0];
  const double num = mi - emi;
  
  if (std::abs(num) < eps) {
    res[0] = 0.0;
  } else {
    const double denom = mn - emi;
    res[0] = std::abs(denom) < eps ? NA_REAL : num / denom;
  }
  
  return res;
}
