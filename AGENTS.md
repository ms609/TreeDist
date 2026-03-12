# TreeDist — Agent Memory

## Project Overview

**TreeDist** is an R package (GPL ≥ 3) providing a suite of topological distance metrics
between phylogenetic trees. The mathematical core is implemented in C++17 and exposed to R
via Rcpp. The primary optimization goal is speed: many real analyses compute pairwise
distances over hundreds or thousands of trees, so inner loops must be tight.

Current version: `2.12.0.9000` (development).  
CRAN package page: <https://cran.r-project.org/package=TreeDist>

---

## Repository Layout

```
TreeDist/
├── src/                  # C++17 source (the main optimization target)
├── R/                    # R wrapper layer and pure-R helpers
├── benchmark/            # Microbenchmark scripts (bench package)
├── tests/testthat/       # Unit tests
├── data-raw/             # Scripts that regenerate lookup tables / data
├── vignettes/            # User-facing tutorials
└── inst/                 # Installed extras
```

**Sibling repository** — `../TreeTools` is a companion package that TreeDist links against
at the C++ level (`LinkingTo: TreeTools`). Edits to TreeTools headers (especially
`SplitList.h`) can affect TreeDist performance and can be pushed to CRAN independently
when ready.

---

## C++ Source Files

| File | Size | Purpose |
|------|------|---------|
| `tree_distances.cpp` | 22 KB | Main distance calculations; calls into CostMatrix / LAP |
| `tree_distances.h` | 15 KB | **CostMatrix** class; cache-aligned storage; `findRowSubmin` hot path |
| `lap.cpp` | 10 KB | Jonker-Volgenant LAPJV linear assignment; extensively hand-optimized |
| `spr.cpp` | 7 KB | SPR distance approximation |
| `spr_lookup.cpp` | — | SPR lookup-table implementation |
| `nni_distance.cpp` | 16 KB | NNI distance approximations; HybridBuffer allocation |
| `li_diameters.h` | 30 KB | Precomputed NNI diameter lookup tables |
| `information.h` | 6 KB | log₂ / factorial lookup tables (max 8192); cached at startup |
| `binary_entropy_counts.cpp` | — | Binary entropy calculations |
| `day_1985.cpp` | 10 KB | Consensus tree information |
| `hmi.cpp` | 6 KB | Hierarchical Mutual Information |
| `hpart.cpp` | 7 KB | Hierarchical partition structures |
| `reduce_tree.cpp` | 11 KB | Prune trees to common tip set before distance calculation |
| `path_vector.cpp` | 3 KB | Path distance vector |
| `mast.cpp` | 5 KB | Maximum Agreement Subtree |
| `RcppExports.cpp` | 20 KB | Auto-generated Rcpp glue (do not edit by hand) |
| `ints.h` | — | Fixed-width integer typedefs (`splitbit`, `int16`, `int32`, …) |

---

## Key Optimizations Already in Place

Understanding what has already been done avoids duplicating effort.

### Memory / cache
- **64-byte alignment** (`alignas(64)`) on `CostMatrix::data_` and `t_data_`.
- `CostMatrix` pads row width to a multiple of `BLOCK_SIZE = 16` so SIMD loads never
  straddle cache lines.
- A **transposed copy** of the cost matrix is maintained to allow column-wise access with
  sequential memory reads.

### Arithmetic
- **Lookup tables** in `information.h` for `log2` (values 0–(SL_MAX_TIPS-1)²) and
  `log2_factorial` (up to 8192); initialized via `__attribute__((constructor))`.
  Using these avoids runtime `std::log2()` calls on hot paths.
- **Loop unrolling** — `CostMatrix::findRowSubmin()` manually unrolls 4× (comment:
  "gives ~20% speedup").
- **`__restrict__`** pointer annotations in `lap.cpp` and `tree_distances.h` to enable
  compiler alias analysis.
- **`__builtin_assume_aligned(ptr, 64)`** hints around inner loops.
- **`__builtin_popcount()`** for split-bit counting.

### Algorithm
- **Column reduction in reverse order** in LAPJV (faster convergence).
- **`findRowSubmin` two-pass strategy** avoids redundant comparisons.
- **LAPJV v2.10.0** achieved a 2× speedup for large matrices via reorganisation.
- **`HybridBuffer<T, StackSize>`** in `nni_distance.cpp`: small allocations go on the
  stack (thresholds: 512 splits, 16 bins) to avoid heap overhead.
- **Tree reduction** (`reduce_tree.cpp`): trees are pruned to their common tip set before
  any distance is computed; this is a major algorithmic win for partially-overlapping
  taxon sets.

### Types
- `cost = int_fast64_t` for LAP to avoid overflow and exploit 64-bit registers.
- `BIG = numeric_limits<cost>::max() / SL_MAX_SPLITS` — the divide avoids integer
  overflow inside the LAP inner loop.
- `ROUND_PRECISION = 2048²` for safe rounding in cost scaling.

---

## Benchmark Infrastructure

All benchmarks live in `benchmark/` and use the `bench` package.

```r
source("benchmark/_init.R")   # loads TreeTools, TreeDist; defines Benchmark()
```

`Benchmark(expr)` wraps `bench::mark()` with `min_iterations = 3` and `time_unit = "us"`.
When run non-interactively (e.g. CI) results are serialised to `.bench.Rds` files.
`_compare_results.R` compares PR results against `main` and fails the build on a
regression > 4 %.

### ⚠ Benchmarking protocol — always use a release install

**Do not benchmark with `devtools::load_all()`.**  It appends
`-UNDEBUG -Wall -pedantic -g -O0` to the compiler flags, which overrides the
`-O2` in `~/.R/Makevars` and produces an unrepresentative (typically 3–5×
slower) build.  All timing figures in this file were produced from a release
install unless noted otherwise.

**Also beware of stale `.o` files.**  `devtools::load_all()` writes debug-
compiled object files to `src/`.  A subsequent `install.packages()` sees
up-to-date timestamps and skips recompilation, silently producing a "release"
install that contains `-O0` objects.  `bench_ab()` and `bench_release()` now
clean `src/*.o` and `src/*.dll` before installing to prevent this.

#### Preferred workflow: A/B comparison (`compare-ab.R`)

`compare-ab.R` installs a renamed copy of the package as **TreeDistRef**,
then installs the current working-tree source as **TreeDist**, and loads
both into a single `Rscript` subprocess for `bench::mark()` comparison.
Because both packages share the same process, CPU state / thermals /
background load affect each measurement equally — eliminating the
environmental noise that plagued the old `bench_release()` approach.

```r
source("benchmark/compare-ab.R")

# 1. On the code you want as baseline:
install_ref()          # installs as TreeDistRef to benchmark/ref_lib/

# 2. Make your changes, then:
bench_ab()             # installs dev TreeDist, compares both in subprocess
```

`install_ref()` copies the source, renames the package (DESCRIPTION,
NAMESPACE, RcppExports symbol prefixes, `R_init_` entry point), cleans
stale `.o` files, and installs to `benchmark/ref_lib/`.  `bench_ab()`
does the same for the dev version (to a temp lib), then runs the
comparison.  Results are saved to `benchmark/results/ab.Rds`.

Note: `compare-ab.R` and `compare-release.R` are local comparison tools
and intentionally do **not** start with `bench-` — the `bench-` prefix
is reserved for benchmark test scripts run by CI (`_run_benchmarks.R`).

MSD (MatchingSplitDistance) serves as a built-in **canary**: it does not
use `add_ic_element`, so any difference on MSD indicates a build-level
or environmental problem rather than a real code change.

#### Legacy workflow: `bench_release()` / `bench_compare()`

Still available via `benchmark/compare-release.R`.  Useful for recording
named snapshots, but comparisons across runs are vulnerable to
environmental drift (~15–20% between sessions on this machine).

```r
source("benchmark/compare-release.R")

# 1. Benchmark HEAD (before your patch)
bench_release(label = "baseline")

# 2. Apply your changes, then benchmark again
bench_release(label = "my-optimisation")

# 3. Compare
bench_compare("baseline", "my-optimisation")
```

`bench_release()` installs the current working-tree source to a private temp
library via `install.packages(..., type = "source")` (so Makevars flags apply
correctly) and runs the suite in a `Rscript` subprocess.  Results are saved to
`benchmark/results/<label>.Rds` and persist across sessions.

### Existing benchmark scripts

| Script | What it tests |
|--------|---------------|
| `bench-tree-distances.R` | CID, PID, RF, MCI on 100×50-tip and 40×200-tip tree sets |
| `bench-LAPJV.R` | LAPJV on 40×40, 400×400, 1999×1999 uniform random matrices |
| `bench-PathDist.R` | PathDist on 6×182-tip trees |
| `bench-MCI-openmp.R` | MCI/CID OpenMP scaling on 100×50, 40×200, and 200×200-tip tree sets |

To add a new benchmark, create `benchmark/bench-<topic>.R` following the existing pattern
and add it to `_run_benchmarks.R`.

---

## Profiling with VTune

VTune 2025.9 is installed at:

```
C:\Program Files (x86)\Intel\oneAPI\vtune\2025.9\bin64\vtune.exe
```

Typical workflow:
1. Build TreeDist with debug symbols but optimisations enabled:
   `PKG_CXXFLAGS="-O2 -g"` in `~/.R/Makevars`.
2. Write a driver `.R` script that exercises the hot path in a loop (use the benchmark
   scripts as a starting point).
3. Run a hotspot collection:
   ```
   vtune -collect hotspots -result-dir vtune-out -- Rscript path/to/driver.R
   ```
4. View the report:
   ```
   vtune -report hotspots -result-dir vtune-out
   ```
   Or open the `.vtune` project in the VTune GUI.
5. Pay attention to the **Memory Access** and **Threading** analyses for cache-miss and
   false-sharing diagnostics.

---

## TreeTools Dependency

`../TreeTools` is available locally and editable. It is linked at the C++ level via
`LinkingTo`; the key header consumed by TreeDist is `<TreeTools/SplitList.h>`.

Important constants defined there that affect TreeDist:
- `SL_MAX_TIPS` — maximum number of leaf taxa per tree.
- `SL_MAX_SPLITS` — maximum number of splits.
- `splitbit` — the unsigned integer type used for bitset representation of splits.

If a bottleneck traces back to a TreeTools header (e.g. SplitList layout, bit-width of
`splitbit`), changes can be made in `../TreeTools`, tested locally by re-installing, and
pushed to CRAN when stable.

---

## Development Workflow

```r
# Build and reload (from R)
devtools::load_all()          # fast incremental rebuild
devtools::test()              # run testthat suite

# Or from the shell
R CMD build .
R CMD check TreeDist_*.tar.gz

# Run a single benchmark interactively
source("benchmark/_init.R")
source("benchmark/bench-tree-distances.R")
```

C++ compilation flags are controlled by `src/Makevars.win` (Windows) / `src/Makevars`.
The package requires **C++17** (`CXX_STD = CXX17`).

---

## Completed Optimizations (this dev cycle)

### OpenMP parallelism for pairwise distances (`src/pairwise_distances.cpp`)
Added `#pragma omp parallel for schedule(dynamic)` over the pairwise loop for
`MutualClusteringInfo` / `ClusteringInfoDistance`.  Build infrastructure:

- `src/Makevars` and `src/Makevars.win` created (these did not exist before);
  both set `CXX_STD = CXX17` and `PKG_CXXFLAGS/PKG_LIBS = $(SHLIB_OPENMP_CXXFLAGS)`.
  On platforms without OpenMP the flags are empty and the package builds
  single-threaded.
- `lap()` in `lap.cpp` gained an `allow_interrupt = true` parameter; the batch
  path passes `false` to avoid calling `Rcpp::checkUserInterrupt()` from a
  worker thread.
- `add_ic_element` moved from `tree_distances.cpp` (inside `namespace TreeDist`,
  inaccessible to other TUs) into `tree_distances.h` as a proper `inline`,
  fixing a latent ODR issue and making it visible to `pairwise_distances.cpp`.
- `SplitList` is not move-constructible; `std::vector<SplitList>` therefore
  cannot be used (reserve/emplace_back require movability).  Fix: use
  `std::vector<std::unique_ptr<SplitList>>` instead.
- R fast path: `.SplitDistanceAllPairs()` in `tree_distance_utilities.R` detects
  `Func == MutualClusteringInfoSplits` with a same-tip-set, no-cluster call and
  routes it to `cpp_mutual_clustering_all_pairs()`.

**Benchmark script**: `benchmark/bench-MCI-openmp.R`
**To measure speedup**: install with `install.packages(".", repos=NULL, type="source")`
into a fresh library (avoids `devtools::load_all()` debug flags and Windows DLL lock).
Use `TreeDist:::cpp_mutual_clustering_all_pairs` when calling from an installed package.

#### Measured speedups (release build, -O2, 16-core Windows machine)

| Scenario | Serial R loop | R parallel (16 workers) | OpenMP batch |
|---|---|---|---|
| 100 trees × 50 tips (4 950 pairs) | 117 ms | 78 ms | **17 ms** |
| 40 trees × 200 tips (780 pairs) | 292 ms | 3 350 ms | **35 ms** |

OpenMP vs serial: **7–8×**.  OpenMP vs R-parallel cluster: **5–95×** (R-parallel
incurs ~2–3 s serialisation overhead regardless of problem size).

#### R parallel cluster crossover (50-tip trees, 16 workers)

The R parallel cluster is only competitive with the serial loop when the problem
is large enough to amortise its ~2–3 s IPC/serialisation overhead.  For 50-tip
trees that crossover is around **500 trees (~125 000 pairs)**, where serial takes
~2.9 s and parallel takes ~2.8 s.  The OpenMP batch path remains faster than
both at every size measured.

Practical implication: **`StartParallel()` provides no benefit for MCI/CID** on
machines where OpenMP is available, and actively harms performance for typical
analysis sizes (≤ 200 trees).  The fast path therefore bypasses the R cluster
entirely when `cluster` is `NULL`.

#### Per-call heap allocation elimination via `LapScratch` (attempted, reverted)

Attempted to eliminate the ~8 per-pair `std::vector` / `CostMatrix` heap
allocations in `pairwise_distances.cpp` by defining a `LapScratch` struct that
holds all scratch storage and is reused across pairs (one instance per OpenMP
thread, indexed by `omp_get_thread_num()`).

**Release-build (`-O2`) benchmarks showed no measurable difference (−0.7% to
+0.5% across all scenarios)** — the allocator is fast enough at `-O2` that
allocation overhead is buried in the Dijkstra phase, which dominates throughout.

The implementation was reverted because it added substantial complexity for
zero measured benefit.

#### Fused IC calculation in `mutual_clustering_score` (attempted, reverted)

Attempted to reduce integer arithmetic in the `add_ic_element` inner loop by
replacing the 4 independent `add_ic_element()` calls (8 integer multiplies
total: `nkK * n_tips` + `nk * nK` per call) with a fused inline block
requiring only 2 multiplies.  The remaining 6 products were derived via
addition/subtraction from precomputed `na_n = na * n_tips`, `nA_n`,
`nb_n_arr[]`, `n_sq`, and `den_ab = na * nb`.  Correctness was verified
(max |batch − serial| ≈ 1e-14).

**Release-build (`-O2`) benchmarks showed no reliable improvement (−13% to
+10% across scenarios, within run-to-run noise).**  The bottleneck in the IC
loop is **memory-bound** (random lookups into the `lg2[]` table), not
compute-bound.  Modern x86 integer multiply has 1-cycle throughput, so
saving 6 multiplies per split pair (~6 cycles) is invisible next to lg2
table cache misses (~4 cycles per L1 hit, much more on L2/L3).

Key lesson: `add_ic_element` and the cost-matrix filling loop are dominated
by `lg2[]` table access, not by integer arithmetic.  Future optimisations
should target the table's working set size or access pattern rather than
the surrounding arithmetic.

#### lg2 table shrink: 32 MB → 16 KB via log2 decomposition

The `lg2[]` lookup table stored `log2(i)` for `i = 0..(SL_MAX_TIPS−1)²`,
requiring 4.2 M doubles = **32 MB**.  `add_ic_element` accessed it with
product indices `lg2[nkK * n_tips]` and `lg2[nk * nK]`, so the working set
for 200-tip trees was ~320 KB (indices up to 40 000) — too large for L1.

**Fix:** decompose `log2(a × b) = log2(a) + log2(b)`:

```cpp
// Before: lg2[nkK * n_tips] - lg2[nk * nK]       (indices up to n²)
// After:  lg2[nkK] + lg2_n  - lg2[nk] - lg2[nK]  (indices ≤ n_tips)
```

The table now has only `SL_MAX_TIPS + 1` entries (2 049 doubles = **16 KB**),
trivially fitting in L1 for any tree size.  `lg2_n = lg2[n_tips]` is
precomputed once per distance call and passed to `add_ic_element` as a
new parameter.

**Files changed:**
- `tree_distances.h`: table declaration shrunk; `add_ic_element` gains
  `lg2_n` parameter and uses decomposed lookups.
- `tree_distance_functions.cpp`: `LG2_SIZE` reduced to `SL_MAX_TIPS + 1`.
- `pairwise_distances.cpp`, `tree_distances.cpp`: callers updated to
  precompute and pass `lg2_n`.

**A/B benchmark (release build, same-process comparison via `compare-ab.R`):**

| Scenario | ref | dev | Change |
|---|---|---|---|
| CID 100 × 50-tip | 75.4 ms | 71.2 ms | **−5.6%** |
| CID 40 × 200-tip | 278 ms | 265 ms | **−4.7%** |
| MSD 100 × 50-tip (canary) | 85.0 ms | 85.6 ms | +0.7% |
| MSD 40 × 200-tip (canary) | 402 ms | 407 ms | +1.2% |

Numerical accuracy unchanged (max |ref − dev| ≈ 5.7 × 10⁻¹⁴).

#### OpenMP rollout to all remaining distance metrics (this dev cycle)

Extended OpenMP batch computation to all LAP-based and RF-info metrics.
Added score-only static functions and `cpp_*_all_pairs` exports in
`pairwise_distances.cpp` for:

| C++ batch function | R `Func` intercepted | LAP? |
|---|---|---|
| `cpp_rf_info_all_pairs` | `InfoRobinsonFouldsSplits` | No |
| `cpp_msd_all_pairs` | `MatchingSplitDistanceSplits` | Yes |
| `cpp_msi_all_pairs` | `MatchingSplitInfoSplits` | Yes |
| `cpp_shared_phylo_all_pairs` | `SharedPhylogeneticInfoSplits` | Yes |
| `cpp_jaccard_all_pairs` | `NyeSplitSimilarity`, `JaccardSplitSimilarity` | Yes |

The R fast-path in `.SplitDistanceAllPairs()` was refactored from a single
`if` block into a unified `if (!is.na(nTip) && is.null(cluster))` guard with
per-function branches, returning early only when a batch function matches.
For `JaccardSplitSimilarity`, `k` and `allowConflict` are extracted from `...`
with sensible defaults (1.0, TRUE) if absent.

---

## LAP Profiling Notes

### Scaling behaviour (debug build, uniform random matrices)

| n | median (μs) | implied α |
|---|---|---|
| 10 | 14.8 | — |
| 25 | 28.1 | 0.70 |
| 47 | 76.6 | 1.59 |
| 97 | 377 | 2.20 |
| 197 | 1798 | 2.21 |
| 400 | 14953 | 2.99 |

For the actual tree-distance workload (n ≈ 47 for 50-tip trees, n ≈ 197 for
200-tip trees), the solver is already well into super-quadratic territory.
The Dijkstra augmentation phase is the dominant cost from n ≈ 100 upward.

### Phase analysis

- **Column reduction** (n² total, O(n) per column via transposed sequential scan):
  already well-optimised; the one-time transpose pays for itself immediately.
- **Reduction transfer** (n² total): the `if (j == j1) continue;` branch inside
  the minimum-search loop was investigated as a vectorisation blocker. Assembly
  analysis (GCC 14 / MinGW, `-O2 -march=native`) showed that **neither** the
  branched nor the split form vectorises — `int_fast64_t` has no SIMD vector
  type on this toolchain, and the outer `matches[i] != 1` check creates a
  multi-loop nest the vectoriser refuses regardless. No change was made.
- **Augmenting row reduction** (2 × free_rows × n): calls `findRowSubmin`,
  already 4× manually unrolled.  `nontrivially_less_than()` was called twice
  per free row with identical arguments; **fixed** by caching as `strictly_less`.
- **Augment solution / Dijkstra** (free_rows × n²): the dominant phase for
  n ≥ 100.  The inner update loop iterates over `col_list[up..dim-1]` via
  indirect `j = col_list[k]` indexing, then accesses `row_i[j]`, `v_ptr[j]`,
  and `d[j]` as scattered gathers — preventing SIMD vectorisation.

### Dijkstra restructuring opportunity (unimplemented — needs release-build VTune)

Replace the `col_list` permutation with a `bool scanned[dim]` mask and iterate
directly over `0..dim-1`.  Pros: sequential access to `row_i`, `v_ptr`, `d`
→ auto-vectorisable; no indirect gather.  Cons: visits all `dim` columns each
iteration (vs progressively fewer with `col_list`), so ~2× more comparisons in
the best case.  Net benefit is uncertain and must be measured on a release build.

**To profile:** build with `PKG_CXXFLAGS="-O2 -g -fno-omit-frame-pointer"` and
run `benchmark/vtune-driver.R` through VTune hotspot collection.

### Test suite fix

Added `tests/testthat/setup.R` that opens a null PDF device for the duration of
all test runs.  This suppresses bare `plot()` / `TreeDistPlot()` calls in tests
from appearing in the interactive graphics device, and prevents vdiffr snapshot
rendering from leaking to screen.  vdiffr opens its own `svglite` device on top
of the null device, so snapshot tests are unaffected.

---

## Known Optimization Opportunities / TODOs

- `information.h` `MAX_FACTORIAL_LOOKUP`: **DONE** — increased from 8192 to 32768;
  `lgamma()` fallback retained for values beyond the table (negligible memory cost).
- `information.h` lines 120–122: `log2_factorial_table` is a verbatim copy from
  `TreeSearch/src/expected_mi.cpp`; should be defined once (TreeTools?) and shared.
- LAP inner loop: the 4× manual unroll works well; investigate whether AVX2 SIMD
  intrinsics (`_mm256_*`) could replace it on modern hardware.
- `CostMatrix` transpose: currently maintained as a full second copy; a cache-oblivious
  blocking scheme (already partially implemented via `BLOCK_SIZE`) could reduce memory
  bandwidth further.
- SPR distance (`spr.cpp`, `spr_lookup.cpp`): the algorithm is relatively recent
  (v2.8.0); profiling under VTune may reveal further hot spots.
- OpenMP for other metrics: **DONE** — see "Completed Optimizations" above.
- lg2 table working set: **DONE** — shrunk from 32 MB to 16 KB; see above.
- Large-tree path (`int32` migration, v2.12.0 dev): **AUDITED** — no issues.
  Fix applied: `split_size` vector in `calc_consensus_info()` moved outside the
  Fix applied: `split_size` vector in `calc_consensus_info()` moved outside the
  per-pair loop to avoid repeated heap allocation for large trees.
### Dijkstra restructuring in LAP (attempted, reverted)

Replaced the `col_list` permutation in the Dijkstra augment-solution phase with a
`scanned[]` (`uint8_t`) mask for sequential memory access (enabling hardware prefetch
and potential auto-vectorization of the `row_i[j] - v_ptr[j] - h` loop).

**A/B benchmark result:** Neutral — within ±3% across all scenarios (CID/MSD,
50-tip/200-tip trees). The sequential-access benefit is cancelled by the algorithmic
overhead of visiting all `dim` columns every iteration (vs progressively fewer columns
with `col_list`). Reverted to avoid adding complexity for no measured gain.

Key lesson: the Dijkstra inner loop bottleneck is not memory access pattern but the
number of comparisons performed. The `col_list` scheme's progressive shrinkage of the
unscanned set is the dominant factor at these problem sizes (n ≈ 47–197).

### VTune hotspot analysis (vtune-lap/)

A VTune hotspot collection (`vtune-lap/`) was run on the pre-lg2-shrink build.
Top TreeDist hotspots (approximate CPU time):

| Function | CPU Time | Notes |
|---|---|---|
| `lap` | 0.975s | LAP solver — expected dominant cost |
| `TreeDist::spi_overlap` | 0.671s | **Surprising — larger than lap on this workload** |
| `std::vector<cost>::vector` constructor | 0.623s | CostMatrix allocation per pair |
| `TreeDist::one_overlap` | 0.441s | Called from spi_overlap |
| `TreeDist::add_ic_element` | 0.252s | Since reduced by lg2 shrink |

**Important caveat**: this profile predates the lg2 table shrink and LapScratch
optimizations; relative weights have since shifted.

### spi_overlap / one_overlap optimisation (attempted, reverted)

Based on the VTune profile, attempted two changes to `tree_distances.h`:

1. **`one_overlap` branch simplification**: unified the `a < b` and `a > b` cases
   into a single `lo = min(a,b)` / `hi = max(a,b)` path (retains the `a == b`
   special case). This removes one unpredictable branch. **Kept** — cleaner code,
   neutral performance.

2. **`spi_overlap` single-pass accumulation**: replaced the original 4-pass loop
   structure (3 while-loops with pointer resets + 1 for-loop) with a single pass
   that accumulates all four bitwise conditions (`a&b`, `~a&b`, `a&~b`, `~(a|b)`)
   at once, applying the last-bin mask only to `~(a|b)`. Also tried adding an
   early-exit `if (all four set) return 0` inside the non-last-bin loop.

   **A/B benchmark results (PID = PhylogeneticInfoDistance, release build):**

   | Version | PID 50-tip (ref 112ms) | PID 200-tip (ref 295ms) | CID canary |
   |---|---|---|---|
   | Single-pass, no early exit | 101ms (−10%) | 341ms (+16%) | neutral |
   | Single-pass + early exit | 134ms (+19%) | 431ms (+44%) | neutral |

   Neither version is a net win. The root cause: the per-pass early exits in the
   original code commonly fire at bin 0 for large trees (n_bins=4), so the original
   loads bin 0 four times (register-hot after the first access). The single-pass
   variant always loads ALL n_bins bins, reading more distinct memory locations.

   More importantly, `spi_overlap` cost-matrix filling is a **minority of per-pair
   time** in the current build — LAP dominates (~70–90% at n≈47), so even a 20%
   reduction in spi_overlap would yield only 2–6% overall improvement.

   **Reverted** to the original 4-pass logic. `one_overlap` simplification retained.

Key lesson: the spi_overlap VTune profile was representative of the pre-lg2-shrink
build where the add_ic_element inner loop was slower.  The current LAP is the
overwhelmingly dominant bottleneck for PID and CID.
