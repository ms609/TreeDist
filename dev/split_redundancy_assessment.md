# Global unique-split index (issue #57): assessment

Question: are the two outstanding items in #57 worth building?

1. A global index of unique splits across all trees (k-way merge), each tree
   stored as links into it.
2. A per-unique-split property cache, and a cache of comparisons between
   pairs of unique splits.

**Conclusion: no for MCI / CID, the headline metric. Close the remaining
items as won't-fix.** The one place with real headroom is
`SharedPhylogeneticInfo()` on redundant tree sets, but there the larger cost
is the full LAP, which a split index does not touch.

## Method

- `dev/split_redundancy.R` builds tree sets of decreasing redundancy
  (independent random trees; an NNI random walk with 10 or 2 moves per step,
  like a thinned MCMC chain; independent 5-NNI perturbations of one tree)
  at 50 tips × 200 trees, 200 × 100 and 1000 × 30.
- It reports the unique-split ratio, the mean fraction of splits a pair shares
  exactly, the R-side `as.Splits()` time, and the single-threaded run time of
  `cpp_mutual_clustering_all_pairs()` and `cpp_shared_phylo_all_pairs()`.
- `dev/split_redundancy_phase_timers.patch` adds timers around the
  exact-match scan, the cost-matrix fill and the LAP, exposed as
  `cpp_phase_times()` (apply, run `Rcpp::compileAttributes()`, install to a
  separate library, pass that library to the script).
- `dev/split_redundancy_lookup.cpp` times filling a cost matrix by computing
  each MCI cell versus gathering it from a precomputed U × U table.

Machine: 4-core Xeon 2.8 GHz, R 4.3.3, g++ -O2, one thread.

## Results

Phase columns are shares of the C++ call; MCI shares sum to < 1 because
pairs that match completely exit before the fill.

| scenario | tips | trees | unique / total | pair shared | `as.Splits` ms | MCI ms | match | fill | LAP | PID ms | fill | LAP |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| random | 50 | 200 | 0.836 | 0.003 | 6 | 1822 | 0.02 | 0.35 | 0.62 | 1530 | 0.24 | 0.75 |
| walk-10NNI | 50 | 200 | 0.140 | 0.059 | 5 | 1789 | 0.01 | 0.31 | 0.67 | 1851 | 0.24 | 0.76 |
| walk-2NNI | 50 | 200 | 0.037 | 0.253 | 7 | 1213 | 0.02 | 0.35 | 0.62 | 1960 | 0.26 | 0.73 |
| perturb-5NNI | 50 | 200 | 0.022 | 0.813 | 6 | 107 | 0.16 | 0.40 | 0.36 | 821 | 0.42 | 0.57 |
| random | 200 | 100 | 0.969 | 0.001 | 5 | 12873 | 0.00 | 0.22 | 0.78 | 13798 | 0.14 | 0.86 |
| walk-10NNI | 200 | 100 | 0.050 | 0.334 | 7 | 5531 | 0.01 | 0.29 | 0.70 | 9115 | 0.21 | 0.79 |
| walk-2NNI | 200 | 100 | 0.019 | 0.749 | 4 | 607 | 0.03 | 0.38 | 0.57 | 5202 | 0.35 | 0.65 |
| perturb-5NNI | 200 | 100 | 0.025 | 0.951 | 4 | 62 | 0.26 | 0.36 | 0.20 | 3087 | 0.53 | 0.47 |
| random | 1000 | 30 | 0.998 | 0.000 | 13 | 103800 | 0.00 | 0.10 | 0.90 | 129415 | 0.06 | 0.94 |
| walk-10NNI | 1000 | 30 | 0.043 | 0.904 | 11 | 337 | 0.04 | 0.41 | 0.51 | 21945 | 0.33 | 0.67 |
| walk-2NNI | 1000 | 30 | 0.035 | 0.980 | 11 | 33 | 0.22 | 0.22 | 0.12 | 14488 | 0.51 | 0.49 |
| perturb-5NNI | 1000 | 30 | 0.038 | 0.990 | 14 | 26 | 0.29 | 0.12 | 0.02 | 12584 | 0.57 | 0.43 |

Cost-matrix fill, ns per cell (`split_redundancy_lookup.cpp`):

| tips | U | table | compute | lookup |
|---|---|---|---|---|
| 50 | 500 | 2 MB | 11.1 | 1.9 |
| 200 | 1000 | 8 MB | 20.7 | 6.1 |
| 200 | 5000 | 200 MB | 21.9 | 12.8 |
| 1000 | 1500 | 18 MB | 56.7 | 3.6 |
| 1000 | 5000 | 200 MB | 62.4 | 8.6 |

## Reading

- **Per-tree work is already negligible.** `as.Splits()` plus SplitList
  construction costs ≤ 15 ms against seconds of pairwise work, so item 1's
  "load once" goal is met; storing trees as links into a global index would
  save memory, not time.
- **Random / diverse sets: nothing to cache.** The unique ratio is ≈ 1, so
  every pair of unique splits is compared about once; a pair cache would be
  built and read once, and would need (N·n)² entries (≈ 3 GB at 100 × 200
  tips). The LAP is 62–94 % of run time here.
- **Redundant sets, MCI: the exact-match pre-scan already captures it.**
  MCI falls from 104 s (random) to 26–337 ms (redundant) at 1000 tips. What
  remains is fill on the few unmatched splits; a look-up table (3–15× faster
  per cell) would save at most ~0.3–0.4 of an already small time, and the
  exact-match scan (≤ 0.29) would become an integer comparison rather than
  vanish.
- **Redundant sets, PID: some headroom, but the LAP dominates.** PID has no
  exact-match pre-scan (identical splits need not be matched to each other in
  the optimal assignment), so it fills and solves the full matrix even for
  near-identical trees: 12.6 s for 30 near-identical 1000-tip trees versus
  26 ms for MCI. A unique-pair cache would cut the fill share (0.33–0.57) by
  70–90 %, a 1.3–2× bound for PID only. A better target for that case is the
  LAP itself (e.g. warm-starting from the identity assignment on shared
  splits), which is outside #57.

## Follow-up: pinning identical splits in PID

`shared_phylo_pinned()` in `src/pairwise_distances.cpp` matches identical
splits to each other, solves the LAP on the rest, and then proves the
combined assignment optimal with LP duals. Pairs it cannot certify are
released into the residual LAP; the full LAP is the fallback. Scores are
bit-identical to the full LAP.

The trap is real: forcing identical matches without a certificate lowered
PID in 104 of 300 near-identical pairs (8–30 tips, 1–6 NNI moves).

Cells between two pinned splits are never computed. Both splits come from
one tree, so they are compatible, and their cost is bounded below by a
function of the two split sizes (nested or disjoint, whichever is cheaper).
Separately, for distinct compatible splits A and B,
SPI(A, B) <= min(h(A), h(B)) - log2(3), checked for every size pair up to
2048 tips.

Pinning runs when at least 70% of splits are shared. Single thread, min of
two runs, `cpp_shared_phylo_all_pairs()`; run-to-run noise on this machine
is about ±10%:

| scenario | tips × trees | before (s) | after (s) |
|---|---|---|---|
| random | 50 × 100 | 0.40 | 0.37 |
| walk-10NNI | 50 × 100 | 0.39 | 0.44 |
| walk-2NNI | 50 × 100 | 0.27 | 0.28 |
| perturb-5NNI | 50 × 100 | 0.17 | 0.21 |
| random | 200 × 50 | 2.95 | 3.37 |
| walk-10NNI | 200 × 50 | 1.54 | 1.59 |
| walk-2NNI | 200 × 50 | 0.91 | 0.73 |
| perturb-5NNI | 200 × 50 | 0.73 | 0.32 |
| random | 1000 × 15 | 31.8 | 38.1 |
| walk-10NNI | 1000 × 15 | 4.50 | 1.86 |
| walk-2NNI | 1000 × 15 | 3.50 | 0.55 |
| perturb-5NNI | 1000 × 15 | 3.18 | 0.53 |

Random sets never pin, so they run the same code as before. Three
interleaved A/B runs on random sets gave 1.90–2.10 s (before) against
1.97–2.15 s (after) at 200 × 40, and 8.40–8.51 s against 8.20–9.19 s at
1000 × 8, so the gaps in the random rows of the table are mostly noise.

`PhylogeneticInfoDistance()` on the CI benchmark sets takes 0.064 s instead
of 0.130 s (`as.phylo(0:99, 50)`), 0.089 s instead of 0.282 s
(`as.phylo(0:39, 180)`) and 0.311 s instead of 1.146 s
(`as.phylo(0:49, 280)`).
