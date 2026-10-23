#ifndef _TREEDIST_NNI_EXACT_H
#define _TREEDIST_NNI_EXACT_H

#include <Rcpp/Lightest>
#include <cstddef>
#include <cstdint>
#include <vector>

// Largest tree that nni_exact_small() can represent: each split is held as a
// 16-bit mask, and a tree's (n - 3) splits are packed into four 64-bit words.
constexpr int32_t NNI_EXACT_MAX_TIPS = 16;

// Trees visited (from both ends) before a region search gives up and falls
// back to the diameter bound.
constexpr std::size_t NNI_EXACT_BUDGET = std::size_t(1) << 21;

// Exact NNI distance between two binary unrooted trees on `n_tip` leaves,
// each given as its n_tip - 3 non-trivial splits: masks over leaves
// 1..n_tip - 1 (i.e. polarized so that leaf 0 is never set).
// Returns -1 if more than `budget` trees would need to be visited.
int32_t nni_exact_small(const std::vector<uint16_t>& a,
                        const std::vector<uint16_t>& b,
                        int32_t n_tip, std::size_t budget);

// Upper bound on NNI distance: sum over regions of unmatched edges of the
// exact distance between the region's two subtrees, for regions with at most
// `max_tips` leaves; other regions contribute nni_region_bound().
// Edges in postorder.
int32_t nni_region_upper(const Rcpp::IntegerMatrix& edge1,
                         const Rcpp::IntegerMatrix& edge2,
                         int32_t n_tip, int32_t max_tips);

// Contribution of a region with `subtree_edges` unmatched edges to
// best_upper; defined in nni_distance.cpp.
int32_t nni_region_bound(int32_t subtree_edges);

#endif
