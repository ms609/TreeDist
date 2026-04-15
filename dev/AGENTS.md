# TreeDist — Agent Memory

## Project Overview

**TreeDist** is an R package (GPL ≥ 3) providing a suite of topological
distance metrics between phylogenetic trees. The mathematical core is
implemented in C++17 and exposed to R via Rcpp. The primary optimization
goal is speed: many real analyses compute pairwise distances over
hundreds or thousands of trees, so inner loops must be tight.

Current version: `2.12.0.9000` (development).  
CRAN package page: <https://cran.r-project.org/package=TreeDist>

------------------------------------------------------------------------

## Repository Layout

    TreeDist/
    ├── src/                  # C++17 source (the main optimization target)
    ├── R/                    # R wrapper layer and pure-R helpers
    ├── benchmark/            # Microbenchmark scripts (bench package)
    ├── tests/testthat/       # Unit tests
    ├── data-raw/             # Scripts that regenerate lookup tables / data
    ├── vignettes/            # User-facing tutorials
    └── inst/                 # Installed extras

**Sibling repository** — `../TreeTools` is a companion package that
TreeDist links against at the C++ level (`LinkingTo: TreeTools`). Edits
to TreeTools headers (especially `SplitList.h`) can affect TreeDist
performance and can be pushed to CRAN independently when ready.

------------------------------------------------------------------------

## C++ Source Files

| File                        | Size  | Purpose                                                                                                   |
|-----------------------------|-------|-----------------------------------------------------------------------------------------------------------|
| `tree_distances.cpp`        | 22 KB | Main distance calculations; calls into CostMatrix / LAP                                                   |
| `tree_distances.h`          | 6 KB  | Tree-distance scoring functions (`add_ic_element`, `one_overlap`, `spi_overlap`); includes `lap.h`        |
| `lap.h`                     | 15 KB | **CostMatrix** class; `LapScratch`; `lap()` declarations; cache-aligned storage; `findRowSubmin` hot path |
| `lap.cpp`                   | 10 KB | Jonker-Volgenant LAPJV linear assignment; extensively hand-optimized; includes only `lap.h`               |
| `spr.cpp`                   | 7 KB  | SPR distance approximation                                                                                |
| `spr_lookup.cpp`            | —     | SPR lookup-table implementation                                                                           |
| `nni_distance.cpp`          | 16 KB | NNI distance approximations; HybridBuffer allocation                                                      |
| `li_diameters.h`            | 30 KB | Precomputed NNI diameter lookup tables                                                                    |
| `information.h`             | 6 KB  | log₂ / factorial lookup tables (max 8192); cached at startup                                              |
| `binary_entropy_counts.cpp` | —     | Binary entropy calculations                                                                               |
| `day_1985.cpp`              | 10 KB | Consensus tree information                                                                                |
| `hmi.cpp`                   | 6 KB  | Hierarchical Mutual Information                                                                           |
| `hpart.cpp`                 | 7 KB  | Hierarchical partition structures                                                                         |
| `reduce_tree.cpp`           | 11 KB | Prune trees to common tip set before distance calculation                                                 |
| `path_vector.cpp`           | 3 KB  | Path distance vector                                                                                      |
| `mast.cpp`                  | 5 KB  | Maximum Agreement Subtree                                                                                 |
| `RcppExports.cpp`           | 20 KB | Auto-generated Rcpp glue (do not edit by hand)                                                            |
| `ints.h`                    | —     | Fixed-width integer typedefs (`splitbit`, `int16`, `int32`, …)                                            |

------------------------------------------------------------------------

## Key Optimizations Already in Place

Understanding what has already been done avoids duplicating effort.

### Memory / cache

- **64-byte alignment** (`alignas(64)`) on `CostMatrix::data_` and
  `t_data_`.
- `CostMatrix` pads row width to a multiple of `BLOCK_SIZE = 16` so SIMD
  loads never straddle cache lines.
- A **transposed copy** of the cost matrix is maintained to allow
  column-wise access with sequential memory reads.

### Arithmetic

- **Lookup tables** in `information.h` for `log2` (values
  0–(SL_MAX_TIPS-1)²) and `log2_factorial` (up to 8192); initialized via
  `__attribute__((constructor))`. Using these avoids runtime
  `std::log2()` calls on hot paths.
- **Loop unrolling** — `CostMatrix::findRowSubmin()` manually unrolls 4×
  (comment: “gives ~20% speedup”).
- **`__restrict__`** pointer annotations in `lap.cpp` and
  `tree_distances.h` to enable compiler alias analysis.
- **`__builtin_assume_aligned(ptr, 64)`** hints around inner loops.
- **`__builtin_popcount()`** for split-bit counting.

### Algorithm

- **Column reduction in reverse order** in LAPJV (faster convergence).
- **`findRowSubmin` two-pass strategy** avoids redundant comparisons.
- **LAPJV v2.10.0** achieved a 2× speedup for large matrices via
  reorganisation.
- **`HybridBuffer<T, StackSize>`** in `nni_distance.cpp`: small
  allocations go on the stack (thresholds: 512 splits, 16 bins) to avoid
  heap overhead.
- **Tree reduction** (`reduce_tree.cpp`): trees are pruned to their
  common tip set before any distance is computed; this is a major
  algorithmic win for partially-overlapping taxon sets.

### Types

- `cost = int_fast64_t` for LAP to avoid overflow and exploit 64-bit
  registers.
- `BIG = numeric_limits<cost>::max() / SL_MAX_SPLITS` — the divide
  avoids integer overflow inside the LAP inner loop.
- `ROUND_PRECISION = 2048²` for safe rounding in cost scaling.

------------------------------------------------------------------------

## Benchmark Infrastructure

All benchmarks live in `benchmark/` and use the `bench` package.

``` r
source("benchmark/_init.R")   # loads TreeTools, TreeDist; defines Benchmark()
```

`Benchmark(expr)` wraps `bench::mark()` with `min_iterations = 3` and
`time_unit = "us"`. When run non-interactively (e.g. CI) results are
serialised to `.bench.Rds` files. `_compare_results.R` compares PR
results against `main` and fails the build on a regression \> 4 %.

### ⚠ Benchmarking protocol — always use a release install

**Do not benchmark with `devtools::load_all()`.** It appends
`-UNDEBUG -Wall -pedantic -g -O0` to the compiler flags, which overrides
the `-O2` in `~/.R/Makevars` and produces an unrepresentative (typically
3–5× slower) build. All timing figures in this file were produced from a
release install unless noted otherwise.

**Also beware of stale `.o` files.** `devtools::load_all()` writes
debug- compiled object files to `src/`. A subsequent
[`install.packages()`](https://rdrr.io/r/utils/install.packages.html)
sees up-to-date timestamps and skips recompilation, silently producing a
“release” install that contains `-O0` objects. `bench_ab()` and
`bench_release()` now clean `src/*.o` and `src/*.dll` before installing
to prevent this.

#### Preferred workflow: A/B comparison (`compare-ab.R`)

`compare-ab.R` installs a renamed copy of the package as
**TreeDistRef**, then installs the current working-tree source as
**TreeDist**, and loads both into a single `Rscript` subprocess for
`bench::mark()` comparison. Because both packages share the same
process, CPU state / thermals / background load affect each measurement
equally — eliminating the environmental noise that plagued the old
`bench_release()` approach.

``` r
source("benchmark/compare-ab.R")

# 1. On the code you want as baseline:
install_ref()          # installs as TreeDistRef to benchmark/ref_lib/

# 2. Make your changes, then:
bench_ab()             # installs dev TreeDist, compares both in subprocess
```

`install_ref()` copies the source, renames the package (DESCRIPTION,
NAMESPACE, RcppExports symbol prefixes, `R_init_` entry point), cleans
stale `.o` files, and installs to `benchmark/ref_lib/`. `bench_ab()`
does the same for the dev version (to a temp lib), then runs the
comparison. Results are saved to `benchmark/results/ab.Rds`.

Note: `compare-ab.R` and `compare-release.R` are local comparison tools
and intentionally do **not** start with `bench-` — the `bench-` prefix
is reserved for benchmark test scripts run by CI (`_run_benchmarks.R`).

MSD (MatchingSplitDistance) serves as a built-in **canary**: it does not
use `add_ic_element`, so any difference on MSD indicates a build-level
or environmental problem rather than a real code change.

#### Legacy workflow: `bench_release()` / `bench_compare()`

Still available via `benchmark/compare-release.R`. Useful for recording
named snapshots, but comparisons across runs are vulnerable to
environmental drift (~15–20% between sessions on this machine).

``` r
source("benchmark/compare-release.R")

# 1. Benchmark HEAD (before your patch)
bench_release(label = "baseline")

# 2. Apply your changes, then benchmark again
bench_release(label = "my-optimisation")

# 3. Compare
bench_compare("baseline", "my-optimisation")
```

`bench_release()` installs the current working-tree source to a private
temp library via `install.packages(..., type = "source")` (so Makevars
flags apply correctly) and runs the suite in a `Rscript` subprocess.
Results are saved to `benchmark/results/<label>.Rds` and persist across
sessions.

### Existing benchmark scripts

| Script                   | What it tests                                                       |
|--------------------------|---------------------------------------------------------------------|
| `bench-tree-distances.R` | CID, PID, RF, MCI on 100×50-tip and 40×200-tip tree sets            |
| `bench-LAPJV.R`          | LAPJV on 40×40, 400×400, 1999×1999 uniform random matrices          |
| `bench-PathDist.R`       | PathDist on 6×182-tip trees                                         |
| `bench-MCI-openmp.R`     | MCI/CID OpenMP scaling on 100×50, 40×200, and 200×200-tip tree sets |

To add a new benchmark, create `benchmark/bench-<topic>.R` following the
existing pattern and add it to `_run_benchmarks.R`.

------------------------------------------------------------------------

## Profiling with VTune

VTune 2025.9 is installed at:

    C:\Program Files (x86)\Intel\oneAPI\vtune\2025.9\bin64\vtune.exe

Typical workflow: 1. Build TreeDist with debug symbols but optimisations
enabled: `PKG_CXXFLAGS="-O2 -g"` in `~/.R/Makevars`. 2. Write a driver
`.R` script that exercises the hot path in a loop (use the benchmark
scripts as a starting point). 3. Run a hotspot collection:
`vtune -collect hotspots -result-dir vtune-out -- Rscript path/to/driver.R`
4. View the report: `vtune -report hotspots -result-dir vtune-out` Or
open the `.vtune` project in the VTune GUI. 5. Pay attention to the
**Memory Access** and **Threading** analyses for cache-miss and
false-sharing diagnostics.

------------------------------------------------------------------------

## TreeTools Dependency

`../TreeTools` is available locally and editable. It is linked at the
C++ level via `LinkingTo`; the key header consumed by TreeDist is
`<TreeTools/SplitList.h>`.

Important constants defined there that affect TreeDist: - `SL_MAX_TIPS`
— maximum number of leaf taxa per tree. - `SL_MAX_SPLITS` — maximum
number of splits. - `splitbit` — the unsigned integer type used for
bitset representation of splits.

If a bottleneck traces back to a TreeTools header (e.g. SplitList
layout, bit-width of `splitbit`), changes can be made in `../TreeTools`,
tested locally by re-installing, and pushed to CRAN when stable.

------------------------------------------------------------------------

## Development Workflow

``` r
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

C++ compilation flags are controlled by `src/Makevars.win` (Windows) /
`src/Makevars`. The package requires **C++17** (`CXX_STD = CXX17`).

### Documentation and spelling checks

After any work that adds or modifies roxygen comments, Rd files,
NEWS.md, or vignettes, run:

``` r
devtools::check_man()                # regenerates Rd files and checks for issues
spelling::spell_check_package()      # flags potential misspellings
```

Legitimate technical terms, proper nouns, and code identifiers flagged
by the spell checker should be added to `inst/WORDLIST` (one word per
line, alphabetically sorted). Only fix actual typos in the source.

### Code coverage

Check that existing tests cover all new code. The GHA test suite uses
codecov. To check locally:

``` r
cov <- covr::package_coverage()
covr::report(cov)                        # interactive HTML report
covr::file_coverage(cov, "src/pairwise_distances.cpp")  # per-file summary
```

Aim for full line coverage of new C++ and R code. If a new code path is
not exercised by the existing test suite, add targeted tests in
`tests/testthat/`.

### ⚠ Exploratory / risky R code — use a subprocess

When running experimental R code that may be slow, allocate large
objects, or involve complex loops (e.g., hill-climbing over tree space,
brute-force evaluation of many candidate trees), **run it in a
subprocess** rather than the interactive RStudio session. Long-running
or memory-heavy computations in the main session can freeze or crash
RStudio.

``` r
# Write the experiment to a temp script, then run via Rscript:
writeLines(code, tmp <- tempfile(fileext = ".R"))
system2("Rscript", tmp, stdout = TRUE, stderr = TRUE)

# Or use callr for structured subprocess execution:
callr::r(function() { ... }, timeout = 120)
```

This applies especially to prototype algorithm exploration (e.g., CID
hill-climbing over split space) where per-iteration cost is uncertain.

------------------------------------------------------------------------

## Completed Optimizations (this dev cycle)

### OpenMP parallelism for pairwise distances (`src/pairwise_distances.cpp`)

Added `#pragma omp parallel for schedule(dynamic)` over the pairwise
loop for `MutualClusteringInfo` / `ClusteringInfoDistance`. Build
infrastructure:

- `src/Makevars` and `src/Makevars.win` created (these did not exist
  before); both set `CXX_STD = CXX17` and
  `PKG_CXXFLAGS/PKG_LIBS = $(SHLIB_OPENMP_CXXFLAGS)`. On platforms
  without OpenMP the flags are empty and the package builds
  single-threaded.
- `lap()` in `lap.cpp` gained an `allow_interrupt = true` parameter; the
  batch path passes `false` to avoid calling
  `Rcpp::checkUserInterrupt()` from a worker thread.
- `add_ic_element` moved from `tree_distances.cpp` (inside
  `namespace TreeDist`, inaccessible to other TUs) into
  `tree_distances.h` as a proper `inline`, fixing a latent ODR issue and
  making it visible to `pairwise_distances.cpp`.
- `SplitList` is not move-constructible; `std::vector<SplitList>`
  therefore cannot be used (reserve/emplace_back require movability).
  Fix: use `std::vector<std::unique_ptr<SplitList>>` instead.
- R fast path: `.SplitDistanceAllPairs()` in `tree_distance_utilities.R`
  detects `Func == MutualClusteringInfoSplits` with a same-tip-set,
  no-cluster call and routes it to
  [`cpp_mutual_clustering_all_pairs()`](https://ms609.github.io/TreeDist/dev/reference/cpp_mutual_clustering_all_pairs.md).

**Benchmark script**: `benchmark/bench-MCI-openmp.R` **To measure
speedup**: install with
`install.packages(".", repos=NULL, type="source")` into a fresh library
(avoids `devtools::load_all()` debug flags and Windows DLL lock). Use
`TreeDist:::cpp_mutual_clustering_all_pairs` when calling from an
installed package.

#### Measured speedups (release build, -O2, 16-core Windows machine)

| Scenario                          | Serial R loop | R parallel (16 workers) | OpenMP batch |
|-----------------------------------|---------------|-------------------------|--------------|
| 100 trees × 50 tips (4 950 pairs) | 117 ms        | 78 ms                   | **17 ms**    |
| 40 trees × 200 tips (780 pairs)   | 292 ms        | 3 350 ms                | **35 ms**    |

OpenMP vs serial: **7–8×**. OpenMP vs R-parallel cluster: **5–95×**
(R-parallel incurs ~2–3 s serialisation overhead regardless of problem
size).

#### R parallel cluster crossover (50-tip trees, 16 workers)

The R parallel cluster is only competitive with the serial loop when the
problem is large enough to amortise its ~2–3 s IPC/serialisation
overhead. For 50-tip trees that crossover is around **500 trees (~125
000 pairs)**, where serial takes ~2.9 s and parallel takes ~2.8 s. The
OpenMP batch path remains faster than both at every size measured.

Practical implication:
**[`StartParallel()`](https://ms609.github.io/TreeDist/dev/reference/StartParallel.md)
provides no benefit for MCI/CID** on machines where OpenMP is available,
and actively harms performance for typical analysis sizes (≤ 200 trees).
The fast path therefore bypasses the R cluster entirely when `cluster`
is `NULL`.

#### Per-call heap allocation elimination via `LapScratch` (attempted, reverted)

Attempted to eliminate the ~8 per-pair `std::vector` / `CostMatrix` heap
allocations in `pairwise_distances.cpp` by defining a `LapScratch`
struct that holds all scratch storage and is reused across pairs (one
instance per OpenMP thread, indexed by `omp_get_thread_num()`).

**Release-build (`-O2`) benchmarks showed no measurable difference
(−0.7% to +0.5% across all scenarios)** — the allocator is fast enough
at `-O2` that allocation overhead is buried in the Dijkstra phase, which
dominates throughout.

The implementation was reverted because it added substantial complexity
for zero measured benefit.

#### Fused IC calculation in `mutual_clustering_score` (attempted, reverted)

Attempted to reduce integer arithmetic in the `add_ic_element` inner
loop by replacing the 4 independent `add_ic_element()` calls (8 integer
multiplies total: `nkK * n_tips` + `nk * nK` per call) with a fused
inline block requiring only 2 multiplies. The remaining 6 products were
derived via addition/subtraction from precomputed `na_n = na * n_tips`,
`nA_n`, `nb_n_arr[]`, `n_sq`, and `den_ab = na * nb`. Correctness was
verified (max \|batch − serial\| ≈ 1e-14).

**Release-build (`-O2`) benchmarks showed no reliable improvement (−13%
to +10% across scenarios, within run-to-run noise).** The bottleneck in
the IC loop is **memory-bound** (random lookups into the `lg2[]` table),
not compute-bound. Modern x86 integer multiply has 1-cycle throughput,
so saving 6 multiplies per split pair (~6 cycles) is invisible next to
lg2 table cache misses (~4 cycles per L1 hit, much more on L2/L3).

Key lesson: `add_ic_element` and the cost-matrix filling loop are
dominated by `lg2[]` table access, not by integer arithmetic. Future
optimisations should target the table’s working set size or access
pattern rather than the surrounding arithmetic.

#### lg2 table shrink: 32 MB → 16 KB via log2 decomposition

The `lg2[]` lookup table stored `log2(i)` for `i = 0..(SL_MAX_TIPS−1)²`,
requiring 4.2 M doubles = **32 MB**. `add_ic_element` accessed it with
product indices `lg2[nkK * n_tips]` and `lg2[nk * nK]`, so the working
set for 200-tip trees was ~320 KB (indices up to 40 000) — too large for
L1.

**Fix:** decompose `log2(a × b) = log2(a) + log2(b)`:

``` cpp
// Before: lg2[nkK * n_tips] - lg2[nk * nK]       (indices up to n²)
// After:  lg2[nkK] + lg2_n  - lg2[nk] - lg2[nK]  (indices ≤ n_tips)
```

The table now has only `SL_MAX_TIPS + 1` entries (2 049 doubles = **16
KB**), trivially fitting in L1 for any tree size. `lg2_n = lg2[n_tips]`
is precomputed once per distance call and passed to `add_ic_element` as
a new parameter.

**Files changed:** - `tree_distances.h`: table declaration shrunk;
`add_ic_element` gains `lg2_n` parameter and uses decomposed lookups. -
`tree_distance_functions.cpp`: `LG2_SIZE` reduced to
`SL_MAX_TIPS + 1`. - `pairwise_distances.cpp`, `tree_distances.cpp`:
callers updated to precompute and pass `lg2_n`.

**A/B benchmark (release build, same-process comparison via
`compare-ab.R`):**

| Scenario                  | ref     | dev     | Change    |
|---------------------------|---------|---------|-----------|
| CID 100 × 50-tip          | 75.4 ms | 71.2 ms | **−5.6%** |
| CID 40 × 200-tip          | 278 ms  | 265 ms  | **−4.7%** |
| MSD 100 × 50-tip (canary) | 85.0 ms | 85.6 ms | +0.7%     |
| MSD 40 × 200-tip (canary) | 402 ms  | 407 ms  | +1.2%     |

Numerical accuracy unchanged (max \|ref − dev\| ≈ 5.7 × 10⁻¹⁴).

#### OpenMP rollout to all remaining distance metrics (this dev cycle)

Extended OpenMP batch computation to all LAP-based and RF-info metrics.
Added score-only static functions and `cpp_*_all_pairs` exports in
`pairwise_distances.cpp` for:

| C++ batch function           | R `Func` intercepted                           | LAP? |
|------------------------------|------------------------------------------------|------|
| `cpp_rf_info_all_pairs`      | `InfoRobinsonFouldsSplits`                     | No   |
| `cpp_msd_all_pairs`          | `MatchingSplitDistanceSplits`                  | Yes  |
| `cpp_msi_all_pairs`          | `MatchingSplitInfoSplits`                      | Yes  |
| `cpp_shared_phylo_all_pairs` | `SharedPhylogeneticInfoSplits`                 | Yes  |
| `cpp_jaccard_all_pairs`      | `NyeSplitSimilarity`, `JaccardSplitSimilarity` | Yes  |

The R fast-path in `.SplitDistanceAllPairs()` was refactored from a
single `if` block into a unified `if (!is.na(nTip) && is.null(cluster))`
guard with per-function branches, returning early only when a batch
function matches. For `JaccardSplitSimilarity`, `k` and `allowConflict`
are extracted from `...` with sensible defaults (1.0, TRUE) if absent.

------------------------------------------------------------------------

## LAP Profiling Notes

### Scaling behaviour (debug build, uniform random matrices)

| n   | median (μs) | implied α |
|-----|-------------|-----------|
| 10  | 14.8        | —         |
| 25  | 28.1        | 0.70      |
| 47  | 76.6        | 1.59      |
| 97  | 377         | 2.20      |
| 197 | 1798        | 2.21      |
| 400 | 14953       | 2.99      |

For the actual tree-distance workload (n ≈ 47 for 50-tip trees, n ≈ 197
for 200-tip trees), the solver is already well into super-quadratic
territory. The Dijkstra augmentation phase is the dominant cost from n ≈
100 upward.

### Phase analysis

- **Column reduction** (n² total, O(n) per column via transposed
  sequential scan): already well-optimised; the one-time transpose pays
  for itself immediately.
- **Reduction transfer** (n² total): the `if (j == j1) continue;` branch
  inside the minimum-search loop was investigated as a vectorisation
  blocker. Assembly analysis (GCC 14 / MinGW, `-O2 -march=native`)
  showed that **neither** the branched nor the split form vectorises —
  `int_fast64_t` has no SIMD vector type on this toolchain, and the
  outer `matches[i] != 1` check creates a multi-loop nest the vectoriser
  refuses regardless. No change was made.
- **Augmenting row reduction** (2 × free_rows × n): calls
  `findRowSubmin`, already 4× manually unrolled.
  `nontrivially_less_than()` was called twice per free row with
  identical arguments; **fixed** by caching as `strictly_less`.
- **Augment solution / Dijkstra** (free_rows × n²): the dominant phase
  for n ≥ 100. The inner update loop iterates over `col_list[up..dim-1]`
  via indirect `j = col_list[k]` indexing, then accesses `row_i[j]`,
  `v_ptr[j]`, and `d[j]` as scattered gathers — preventing SIMD
  vectorisation.

### Dijkstra restructuring opportunity (unimplemented — needs release-build VTune)

Replace the `col_list` permutation with a `bool scanned[dim]` mask and
iterate directly over `0..dim-1`. Pros: sequential access to `row_i`,
`v_ptr`, `d` → auto-vectorisable; no indirect gather. Cons: visits all
`dim` columns each iteration (vs progressively fewer with `col_list`),
so ~2× more comparisons in the best case. Net benefit is uncertain and
must be measured on a release build.

**To profile:** build with
`PKG_CXXFLAGS="-O2 -g -fno-omit-frame-pointer"` and run
`benchmark/vtune-driver.R` through VTune hotspot collection.

### Test suite fix

Added `tests/testthat/setup.R` that opens a null PDF device for the
duration of all test runs. This suppresses bare
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) /
[`TreeDistPlot()`](https://ms609.github.io/TreeDist/dev/reference/TreeDistPlot.md)
calls in tests from appearing in the interactive graphics device, and
prevents vdiffr snapshot rendering from leaking to screen. vdiffr opens
its own `svglite` device on top of the null device, so snapshot tests
are unaffected.

------------------------------------------------------------------------

## Known Optimization Opportunities / TODOs

- `information.h` `MAX_FACTORIAL_LOOKUP`: **DONE** — increased from 8192
  to 32768; [`lgamma()`](https://rdrr.io/r/base/Special.html) fallback
  retained for values beyond the table (negligible memory cost).
- `information.h` lines 120–122: `log2_factorial_table` is a verbatim
  copy from `TreeSearch/src/expected_mi.cpp`; should be defined once
  (TreeTools?) and shared.
- LAP inner loop: the 4× manual unroll works well; investigate whether
  AVX2 SIMD intrinsics (`_mm256_*`) could replace it on modern hardware.
- `CostMatrix` transpose: currently maintained as a full second copy; a
  cache-oblivious blocking scheme (already partially implemented via
  `BLOCK_SIZE`) could reduce memory bandwidth further.
- SPR distance (`spr.cpp`, `spr_lookup.cpp`): the algorithm is
  relatively recent (v2.8.0); profiling under VTune may reveal further
  hot spots.
- OpenMP for other metrics: **DONE** — see “Completed Optimizations”
  above.
- lg2 table working set: **DONE** — shrunk from 32 MB to 16 KB; see
  above.
- Large-tree path (`int32` migration, v2.12.0 dev): **AUDITED** — no
  issues. Fix applied: `split_size` vector in `calc_consensus_info()`
  moved outside the per-pair loop to avoid repeated heap allocation for
  large trees.
- Hardware POPCNT (`-mpopcnt`): **DONE** — see “Completed Optimizations”
  below. \### Dijkstra restructuring in LAP (attempted, reverted)

Replaced the `col_list` permutation in the Dijkstra augment-solution
phase with a `scanned[]` (`uint8_t`) mask for sequential memory access
(enabling hardware prefetch and potential auto-vectorization of the
`row_i[j] - v_ptr[j] - h` loop).

**A/B benchmark result:** Neutral — within ±3% across all scenarios
(CID/MSD, 50-tip/200-tip trees). The sequential-access benefit is
cancelled by the algorithmic overhead of visiting all `dim` columns
every iteration (vs progressively fewer columns with `col_list`).
Reverted to avoid adding complexity for no measured gain.

Key lesson: the Dijkstra inner loop bottleneck is not memory access
pattern but the number of comparisons performed. The `col_list` scheme’s
progressive shrinkage of the unscanned set is the dominant factor at
these problem sizes (n ≈ 47–197).

### VTune hotspot analysis (vtune-lap/)

A VTune hotspot collection (`vtune-lap/`) was run on the pre-lg2-shrink
build. Top TreeDist hotspots (approximate CPU time):

| Function                                | CPU Time | Notes                                             |
|-----------------------------------------|----------|---------------------------------------------------|
| `lap`                                   | 0.975s   | LAP solver — expected dominant cost               |
| `TreeDist::spi_overlap`                 | 0.671s   | **Surprising — larger than lap on this workload** |
| `std::vector<cost>::vector` constructor | 0.623s   | CostMatrix allocation per pair                    |
| `TreeDist::one_overlap`                 | 0.441s   | Called from spi_overlap                           |
| `TreeDist::add_ic_element`              | 0.252s   | Since reduced by lg2 shrink                       |

**Important caveat**: this profile predates the lg2 table shrink and
LapScratch optimizations; relative weights have since shifted.

### spi_overlap / one_overlap optimisation (attempted, reverted)

Based on the VTune profile, attempted two changes to `tree_distances.h`:

1.  **`one_overlap` branch simplification**: unified the `a < b` and
    `a > b` cases into a single `lo = min(a,b)` / `hi = max(a,b)` path
    (retains the `a == b` special case). This removes one unpredictable
    branch. **Kept** — cleaner code, neutral performance.

2.  **`spi_overlap` single-pass accumulation**: replaced the original
    4-pass loop structure (3 while-loops with pointer resets + 1
    for-loop) with a single pass that accumulates all four bitwise
    conditions (`a&b`, `~a&b`, `a&~b`, `~(a|b)`) at once, applying the
    last-bin mask only to `~(a|b)`. Also tried adding an early-exit
    `if (all four set) return 0` inside the non-last-bin loop.

    **A/B benchmark results (PID = PhylogeneticInfoDistance, release
    build):**

    | Version                    | PID 50-tip (ref 112ms) | PID 200-tip (ref 295ms) | CID canary |
    |----------------------------|------------------------|-------------------------|------------|
    | Single-pass, no early exit | 101ms (−10%)           | 341ms (+16%)            | neutral    |
    | Single-pass + early exit   | 134ms (+19%)           | 431ms (+44%)            | neutral    |

    Neither version is a net win. The root cause: the per-pass early
    exits in the original code commonly fire at bin 0 for large trees
    (n_bins=4), so the original loads bin 0 four times (register-hot
    after the first access). The single-pass variant always loads ALL
    n_bins bins, reading more distinct memory locations.

    More importantly, `spi_overlap` cost-matrix filling is a **minority
    of per-pair time** in the current build — LAP dominates (~70–90% at
    n≈47), so even a 20% reduction in spi_overlap would yield only 2–6%
    overall improvement.

    **Reverted** to the original 4-pass logic. `one_overlap`
    simplification retained.

Key lesson: the spi_overlap VTune profile was representative of the
pre-lg2-shrink build where the add_ic_element inner loop was slower. The
current LAP is the overwhelmingly dominant bottleneck for PID and CID.

------------------------------------------------------------------------

### VTune hotspot collection workflow (updated)

R’s default `DLLFLAGS` includes `-s` (strip-all), which removes all
symbols from the DLL and makes VTune show addresses like
`func@0x340a9ca80` instead of function names. To produce a profiling
build with full symbols:

1.  Ensure `src/Makevars.win` exists with `-fno-omit-frame-pointer` in
    `PKG_CXXFLAGS`.

2.  Override `DLLFLAGS` via `MAKEFLAGS` when installing to
    `.vtune-lib/`:

    ``` r
    old_mf <- Sys.getenv("MAKEFLAGS")
    Sys.setenv(MAKEFLAGS = "DLLFLAGS=-static-libgcc")
    install.packages(".", lib = ".vtune-lib", repos = NULL, type = "source",
                     INSTALL_opts = "--no-multiarch")
    Sys.setenv(MAKEFLAGS = old_mf)
    ```

3.  Run the VTune hotspot collection:

        vtune -collect hotspots -result-dir vtune-current \
              -- Rscript benchmark/vtune-driver.R

    (`benchmark/vtune-driver.R` exercises CID and PID workloads for ~30
    s each.)

**Do not** add `DLLFLAGS` to `src/Makevars.win` — it has no effect there
because R’s `Makeconf` sets `DLLFLAGS` after the package Makevars is
parsed. The `MAKEFLAGS` environment variable overrides it correctly at
the `make` command level.

After profiling, remove `-fno-omit-frame-pointer` from
`src/Makevars.win`.

### VTune hotspot profile — previous build (vtune-current/, pre-exact-match)

**Workload**: CID 50-tip (100 trees, 30 s) + CID 200-tip (40 trees, 30
s) + PID 50-tip (100 trees, 30 s). Build flags:
`-O2 -fopenmp -fno-omit-frame-pointer -mpopcnt`.

| Function                      | CPU Time | %     | Notes                                            |
|-------------------------------|----------|-------|--------------------------------------------------|
| `mutual_clustering_score`     | 40.9s    | 44.8% | CID/MCI per-pair kernel                          |
| `TreeDist::spi_overlap`       | 9.1s     | 10.0% | Split overlap (PID/MSD path)                     |
| `_popcountdi2`                | 8.8s     | 9.7%  | **Software popcount** — eliminated by `-mpopcnt` |
| `[TreeDist.dll]` (unresolved) | 6.6s     | 7.2%  | Likely inlined callees                           |
| `lap`                         | 5.8s     | 6.4%  | LAP solver                                       |

### VTune hotspot profile — current build (vtune-current-post2/, 2026-03-12)

**Workload**: CID 50/200-tip + PID 50/200-tip + MSD 50/200-tip, ~20 s
each (120 s total). Build flags: `-O2 -fopenmp -fno-omit-frame-pointer`
(inline asm popcnt, no `-mpopcnt`). Trees generated via
`as.phylo(0:N, tips)` — high split overlap (“similar trees” scenario).

| Function                  | CPU Time | %     | Notes                                                                 |
|---------------------------|----------|-------|-----------------------------------------------------------------------|
| `mutual_clustering_score` | 31.2s    | 25.4% | CID/MCI per-pair kernel                                               |
| `msd_score`               | 27.0s    | 22.0% | MSD per-pair kernel (includes LAP)                                    |
| `shared_phylo_score`      | 12.9s    | 10.5% | SPI per-pair kernel (includes LAP)                                    |
| `spi_overlap`             | 12.7s    | 10.3% | Split overlap bitwise ops (called by SPI, Jaccard)                    |
| ucrtbase (memcpy/memset)  | 8.1s     | 6.6%  | CostMatrix / vector init                                              |
| `malloc_base`             | 5.3s     | 4.3%  | Heap allocation                                                       |
| `free_base`               | 5.2s     | 4.2%  | Heap deallocation                                                     |
| R internals               | 8.6s     | 7.0%  | R evaluator, GC, symbol lookup                                        |
| `SplitList::SplitList`    | 0.8s     | 0.7%  | Rcpp → SplitList conversion                                           |
| `lap`                     | 0.7s     | 0.6%  | LAP solver (separate symbol; most LAP time is inlined into \*\_score) |
| `LapScratch::ensure`      | 0.2s     | 0.2%  | Scratch buffer resize                                                 |

**Key observations:** 1. **`_popcountdi2` no longer appears** — the
inline-asm hardware POPCNT fully eliminated it. 2. **`lap` as a separate
symbol is only 0.6%** — but LAP execution time is inlined into
`msd_score` and `shared_phylo_score`, which together account for 32.5%.
The LAP solver dominates both of those functions. 3.
**`mutual_clustering_score` dropped from 44.8% to 25.4%** — partly
because the workload now includes MSD and PID (diluting CID’s share),
and partly from the branchless IC hoisting. 4. **malloc/free together =
8.5%** — this is per-pair CostMatrix and vector allocation. Previous
`LapScratch` pooling attempt showed no gain at `-O2`, but the profile
confirms allocation is a measurable cost. 5. **`spi_overlap` at 10.3%**
is consistent with the previous profile. The 4-pass early-exit structure
was later replaced by a popcount-based single-pass approach (see below).
6. **The “similar trees” workload** means MSD benefits heavily from
exact-match detection (solving reduced LAPs of ~3 instead of ~197). On
random trees, `msd_score` would be much larger relative to
`mutual_clustering_score`.

### Branchless IC hoisting in `mutual_clustering_score` (DONE, kept)

The `add_ic_element()` function was called 4× per (ai, bi) split pair in
the CID cost matrix filling loop. Each call contained two branches
(`nkK && nk && nK` and `numerator != denominator`) plus 4 lg2 table
lookups — 16 lookups and 8 branches per pair total.

**Fix:** algebraically expand the 4 `add_ic_element` calls into a single
branchless expression. Key observations: - `lg2[na]`, `lg2[nA]` are
constant for all `bi` in the inner loop → hoisted per `ai` as
`offset_a = lg2_n - lg2[na]`, `offset_A = lg2_n - lg2[nA]`. - `lg2[nb]`,
`lg2[nB]` are constant for each `bi` → computed once per `bi`. -
`lg2[0] = 0` by convention (0 \* log2(0) = 0), so zero overlap counts
contribute 0 without branching. - The `numerator != denominator` check
(independence case) produces `log2(1) = 0` naturally in floating point,
making the branch unnecessary.

Result: 16 lg2 lookups → 4 per pair; 8 branches → 0; 12 int multiplies →
0.

**Files changed:** `pairwise_distances.cpp` (OpenMP batch path),
`tree_distances.cpp` (serial per-pair path).

**A/B benchmark (release build, same-process comparison):**

| Scenario                  | ref     | dev     | Change        |
|---------------------------|---------|---------|---------------|
| CID 100 × 50-tip          | 66.7 ms | 58.9 ms | **−12%**      |
| CID 40 × 200-tip          | 176 ms  | 154 ms  | **−12.5%**    |
| MSD 100 × 50-tip (canary) | 63.3 ms | 63.9 ms | +0.9% (noise) |
| MSD 40 × 200-tip (canary) | 211 ms  | 216 ms  | +2.4% (noise) |

Numerical accuracy: max \|ref − dev\| ≈ 5.7 × 10⁻¹⁴ (unchanged
precision).

**Post-optimization observation:** CID is now **faster than MSD**
per-pair:

| Metric | 50-tip (μs/pair) | 200-tip (μs/pair) |
|--------|------------------|-------------------|
| CID    | 11.9             | 197               |
| MSD    | 12.9             | 271               |

This revealed that MSD’s cost is dominated by LAP on the full matrix,
not matrix filling.

### MSD exact-match detection in batch path (DONE, kept)

CID already had an exact-match shortcut (detecting
identical/complementary splits and solving a reduced LAP). MSD lacked
this, solving the full n_splits × n_splits LAP every time.

**Root cause investigation:** the `as.phylo(0:N, tips)` benchmark trees
share ~95% of splits (mean RF ≈ 5 for 200-tip trees → ~195 of 197 splits
match). CID’s effective LAP dimension was ~3 while MSD’s was ~197. For
truly random trees (ape::rtree), exact matches ≈ 0.

**Fix:** added exact-match detection to `msd_score()` in
`pairwise_distances.cpp`: - Fused the half_tips complement flip into the
popcount loop (eliminates second row pass) - When XOR popcount = 0
(after flip), marks both splits as matched and breaks - Builds a reduced
cost matrix for remaining unmatched splits - Solves the smaller LAP

**A/B benchmark (release build, same-process comparison):**

| Scenario                         | ref     | dev     | Change      |
|----------------------------------|---------|---------|-------------|
| MSD similar 50-tip (4 950 pairs) | 64.2 ms | 24.5 ms | **−62%**    |
| MSD similar 200-tip (780 pairs)  | 212 ms  | 70.9 ms | **−67%**    |
| MSD random 50-tip                | 106 ms  | 104 ms  | −2% (noise) |
| MSD random 200-tip               | 304 ms  | 299 ms  | −2% (noise) |
| CID similar 50 (canary)          | 57.9 ms | 57.4 ms | −1% (noise) |

Numerically exact (max \|ref − dev\| = 0). No regression for random
trees.

MCMC posteriors and bootstrap replicates typically share most splits —
the “similar trees” scenario is the common real-world case.

Exact-match detection was also applied to the MSI, SPI, and Jaccard
batch paths (all LAP-based) in the same session. The detection logic is
identical (XOR popcount = 0 after complement flip); only the “exact
match contribution” differs per metric. Bugs found and fixed during that
extension: MSI flip inconsistency, SPI consistency cleanup, and Jaccard
`allow_conflict=FALSE` guard (see conversation summary for details).

### CostMatrix pooling in batch path (DONE, kept)

The VTune profile showed 8.5% of CPU time in malloc/free and 6.6% in
ucrtbase memset — caused by per-pair allocation and zero-initialisation
of CostMatrix (two `std::vector<cost>` of `dim8² × 8` bytes each; ~692
KB for 200-tip trees).

**Fix:** pool CostMatrix instances in `LapScratch` so they persist
across pairs.

- Made `CostMatrix` resizable: removed `const` from `dim_`/`dim8_`,
  added default constructor and `resize(new_dim)` that only reallocates
  when the new dimension exceeds the current buffer capacity.
- Fixed `padAfterRow()` and `padTrAfterCol()` to fill only up to
  `dim_ * dim8_` (not `data_.end()`), which would overshoot when the
  pool buffer is oversized from a previous larger dimension.
- Added `score_pool` and `small_pool` CostMatrix fields to `LapScratch`.
- Updated all 5 batch score functions (`mutual_clustering_score`,
  `msd_score`, `msi_score`, `shared_phylo_score`, `jaccard_score`) to
  use the pooled matrices.

**Files changed:** `tree_distances.h` (CostMatrix + LapScratch),
`pairwise_distances.cpp`.

**A/B benchmark (release build, same-process comparison):**

| Scenario            | ref (ms) | dev (ms) | Change   |
|---------------------|----------|----------|----------|
| CID similar 50-tip  | 58.7     | 51.7     | **−12%** |
| CID similar 200-tip | 154.3    | 143.9    | **−7%**  |
| MSD similar 50-tip  | 25.4     | 19.3     | **−24%** |
| MSD similar 200-tip | 79.7     | 66.6     | **−16%** |
| PID similar 50-tip  | 79.0     | 69.0     | **−13%** |
| PID similar 200-tip | 194.6    | 189.1    | **−3%**  |
| CID random 50-tip   | 222.4    | 210.8    | **−5%**  |
| CID random 200-tip  | 773.5    | 760.3    | −2%      |
| MSD random 50-tip   | 107.6    | 97.3     | **−10%** |
| MSD random 200-tip  | 315.9    | 300.5    | **−5%**  |

Numerically exact (max \|ref − dev\| = 0). Gains are largest on similar
trees (where exact-match detection makes per-pair compute cheap, so
allocation overhead is a larger fraction) but also meaningful on random
trees (5–10% for 50-tip).

### Header split: `lap.h` / `tree_distances.h` (DONE, kept)

CI benchmarks for PRs \#178 and \#180 showed a consistent ~10%
regression on the standalone
[`LAPJV()`](https://ms609.github.io/TreeDist/dev/reference/LAPJV.md)
benchmark across all matrix sizes (40, 400, 2000), despite `lap.cpp`
being unchanged between branches.

**Root cause:** `lap.cpp` included `tree_distances.h`, which contains
`namespace TreeDist { ... }` with inline scoring functions
(`one_overlap`, `spi_overlap`, `add_ic_element`). Changes to these
functions — even though LAP never calls them — shifted the instruction
layout of the LAP hot loops (`findRowSubmin`, Dijkstra inner loop)
across alignment boundaries, degrading branch prediction and instruction
fetch throughput.

Key evidence: - PR \#178 (posit-optim-2 → main): the **only**
`tree_distances.h` change was a 5-line refactor of `one_overlap()`, yet
LAPJV regressed 8–12%. CostMatrix and LapScratch were byte-for-byte
identical to main. - The regression was consistent across 3+ independent
CI runs and all matrix sizes. - Tree distance metrics improved 20–45% in
the same PR, confirming it wasn’t an environmental issue.

**Fix:** split `tree_distances.h` into two headers:

| Header             | Content                                                                                                           | Included by                   |
|--------------------|-------------------------------------------------------------------------------------------------------------------|-------------------------------|
| `lap.h`            | Types, constants, CostMatrix, LapScratch, `lap()` declarations                                                    | `lap.cpp`, `tree_distances.h` |
| `tree_distances.h` | `#include "lap.h"` + scoring functions (`namespace TreeDist`) + tree-specific globals (`lg2[]`, `ALL_ONES`, etc.) | all other `.cpp` files        |

`lap.cpp` now includes only `lap.h`. All other translation units
continue to include `tree_distances.h` (which pulls in `lap.h` via
`#include`), so they see identical content.

**Files changed:** `src/lap.h` (new), `src/tree_distances.h` (slimmed),
`src/lap.cpp` (`#include "lap.h"` instead of `"tree_distances.h"`).

**A/B benchmark (release build, same-process comparison, local Windows
machine):**

| Scenario         | ref     | dev     | Change      |
|------------------|---------|---------|-------------|
| CID 100 × 50-tip | 50.3 ms | 50.5 ms | neutral     |
| CID 40 × 200-tip | 143 ms  | 143 ms  | neutral     |
| MSD 100 × 50-tip | 19.5 ms | 19.5 ms | neutral     |
| LAPJV 400×400    | 1.24 ms | 1.24 ms | neutral     |
| LAPJV 1999×1999  | 90.6 ms | 93.1 ms | ~3% (noise) |

On the local machine (Windows, GCC 14, `-O2`) the split is
performance-neutral. The CI regression (Ubuntu, GCC, `-O3`) is
alignment-sensitive and platform-specific; the split eliminates the
mechanism by which future scoring function changes can affect LAP code
generation.

**Maintenance burden:** minimal. The split is along a natural boundary —
`lap.h` contains stable LAP infrastructure that changes rarely;
`tree_distances.h` contains actively-developed scoring functions. No
code duplication.

------------------------------------------------------------------------

### R-level overhead investigation for ClusteringInfoDistance (investigated, not yet optimised)

Profiled the R-level overhead in
[`ClusteringInfoDistance()`](https://ms609.github.io/TreeDist/dev/reference/TreeDistance.md)
(CID) to determine whether further gains are available above the C++
layer.

**Breakdown** (CID 100 × 50-tip, ~67ms total):

| Component                                                                      | Time              | Notes                                   |
|--------------------------------------------------------------------------------|-------------------|-----------------------------------------|
| `as.Splits` (batch path)                                                       | ~2.1 ms           | Necessary — converts trees to SplitList |
| `as.Splits` (entropy path)                                                     | ~3.5 ms           | **DUPLICATE** — same trees re-converted |
| Per-tree R dispatch for `ClusteringEntropy.Splits`                             | ~2.2 ms           | `vapply` over N trees                   |
| Other (dispatch, [`structure()`](https://rdrr.io/r/base/structure.html), etc.) | ~0.5 ms           | General R overhead                      |
| **Total R overhead**                                                           | **~8.2 ms (12%)** |                                         |

**Root cause:**
[`ClusteringInfoDistance()`](https://ms609.github.io/TreeDist/dev/reference/TreeDistance.md)
calls two separate paths: 1.
[`MutualClusteringInfo()`](https://ms609.github.io/TreeDist/dev/reference/TreeDistance.md)
→ `.SplitDistanceAllPairs()` → `as.Splits()` + C++ batch 2.
`.MaxValue(tree1, tree2, ClusteringEntropy)` → `as.Splits()` **again** +
per-tree `ClusteringEntropy.Splits` via `vapply`

The duplicate `as.Splits()` (3.5ms) and per-tree R dispatch (2.2ms)
together account for ~5.7ms (~8.5% of total).

**Potential fixes** (in order of impact/complexity):

1.  **Avoid duplicate `as.Splits()`** — pass pre-computed splits from
    batch path to entropy computation (~5% speedup, simple R-level
    change)
2.  **Add C++ `cpp_clustering_entropy_batch()`** — compute per-tree
    entropy in one C++ call instead of per-tree R dispatch (~3%
    additional speedup)
3.  **Fuse into batch function** — return per-tree entropies as
    attribute alongside pairwise scores (no separate normalisation pass;
    most complex)

**Cost-benefit note:** the overhead is modest in absolute terms (~6ms
for 50-tip, ~4ms for 200-tip) and scales O(N) rather than O(N²). It
matters most for high N (e.g., 1000 × 50-tip trees: ~60ms overhead vs
~500ms compute).

------------------------------------------------------------------------

### R-level fast paths for \*Distance functions (DONE, kept)

Added `.FastDistPath()` helper and `.PairwiseSums()` to
`tree_distance.R`. These pre-compute splits once via `as.Splits()` and
use them for both the C++ batch distance computation AND per-tree
entropy/info computation, avoiding duplicate `as.Splits()` calls.

Applied to 4 distance functions: - `ClusteringInfoDistance`
(`tree_distance_info.R`) — uses `ClusteringEntropy.Splits` -
`DifferentPhylogeneticInfo` (`tree_distance_info.R`) — uses
`SplitwiseInfo.Splits` - `MatchingSplitInfoDistance`
(`tree_distance_msi.R`) — uses `SplitwiseInfo.Splits` -
`InfoRobinsonFoulds` (`tree_distance_rf.R`) — uses
`SplitwiseInfo.Splits`

The fast path fires when: `tree2` is NULL (all-pairs),
`!reportMatching`, same tip set, `nTip >= 4`, no R parallel cluster
configured.

**devtools::load_all() benchmarks (CID, 100 × 50-tip similar trees):** -
Fast path: ~50 ms (was ~60 ms with slow path) → ~17% speedup - For
200-tip trees: ~144 ms (was ~153 ms) → ~6% speedup

### Cross-pairs (ManyMany) batch path (DONE, kept)

Added 6 C++ `cpp_*_cross_pairs()` functions to `pairwise_distances.cpp`:
`mutual_clustering`, `rf_info`, `msd`, `msi`, `shared_phylo`, `jaccard`.
Each accepts two lists of SplitLists and returns an nA × nB matrix.
OpenMP parallelism, CostMatrix pooling, and `LapScratch` reuse all
apply.

R-level batch detection added to `.SplitDistanceManyMany()` in
`tree_distance_utilities.R`, following the same `identical(Func, ...)`
pattern as `.SplitDistanceAllPairs()`.

### SPI / MSI exact-match detection bug (FOUND AND FIXED)

**Bug:** The `shared_phylo_score` and `msi_score` batch functions in
`pairwise_distances.cpp` used greedy exact-match detection (matching
identical splits before solving the LAP), identical to the MCI and MSD
paths. This is **incorrect for SPI and MSI** because `spi_overlap(A, B)`
where B contains A can EXCEED `spi_overlap(A, A)` for the identical
split. The full LAP can therefore find a better global assignment by NOT
matching identical splits.

This produced up to ~1% error (max \|batch − serial\| ≈ 14 on scores ≈
1600). The bug existed since the exact-match detection was added to
SPI/MSI batch paths (earlier in this dev cycle).

MCI and MSD exact-match detection is CORRECT and retained: for MCI,
identical splits provide maximum mutual information; for MSD, identical
splits have distance 0 (provably optimal).

**Fix:** Removed exact-match detection from `shared_phylo_score` and
`msi_score`. These now always solve the full LAP. This is a correctness
fix that costs some performance on similar trees (SPI/MSI now solve the
full n_splits × n_splits LAP instead of a reduced one). MCI, MSD, and
Jaccard batch paths are unaffected.

**Verification:** max \|batch − serial\| = 0 for all metrics across 10 ×
8 cross-pairs and 5 × 5 all-pairs tests.

------------------------------------------------------------------------

### C++ batch entropy/info functions (DONE, kept)

Replaced per-tree `vapply(splits_list, InfoSplitsFn, double(1L))` in
`.FastDistPath()` with two C++ batch functions in
`pairwise_distances.cpp`:

- **`cpp_clustering_entropy_batch(splits_list, n_tip)`**: computes
  per-tree clustering entropy (binary entropy of split proportions, sum
  over splits). Used by `ClusteringInfoDistance`.
- **`cpp_splitwise_info_batch(splits_list, n_tip)`**: computes per-tree
  splitwise information
  (`Log2Unrooted(n) - Log2Rooted(k) - Log2Rooted(n-k)` summed over
  splits). Uses a precomputed `Log2Rooted` table. Used by
  `PhylogeneticInfoDistance`, `MatchingSplitInfoDistance`,
  `InfoRobinsonFoulds`.

**Micro-benchmark (100 trees × 50-tip, debug build):**

| Function          | R vapply | C++ batch | Speedup |
|-------------------|----------|-----------|---------|
| ClusteringEntropy | 2.64 ms  | 0.58 ms   | 4.6×    |
| SplitwiseInfo     | 17.84 ms | 0.40 ms   | **45×** |

SplitwiseInfo was especially slow in R because `SplitwiseInfo.Splits`
calls `vapply(inSplit, Log2Rooted.int, 0)` internally — a nested R
dispatch loop.

**Numerical accuracy:** ClusteringEntropy: max \|diff\| ≈ 1.4e-14
(identical precision). SplitwiseInfo: max \|diff\| ≈ 3.7e-10 for 200-tip
trees (cumulative log2 summation in the C++ table vs TreeTools’
precomputed lookup).

**Files changed:** `src/pairwise_distances.cpp` (new functions),
`R/tree_distance.R` (`.FastDistPath` parameter change),
`R/tree_distance_info.R`, `R/tree_distance_msi.R`,
`R/tree_distance_rf.R` (callers updated).

------------------------------------------------------------------------

### Sort+merge exact-match pre-scan (DONE, kept)

The fused O(n²) exact-match detection + cost-matrix fill loop in CID,
MSD, and Jaccard was the dominant per-pair cost for similar trees. For
200-tip similar trees, ~194 of 197 splits match exactly, so 99% of the
O(n²) loop iterations produced cost-matrix entries that were never used
by the LAP.

**Fix:** two-phase approach:

1.  **Phase 1 — O(n log n) sort+merge matching**: `find_exact_matches()`
    helper canonicalises each split (flip if bit 0 is unset), sorts both
    index arrays by canonical form, then merge-scans to find exact
    matches.
2.  **Phase 2 — O(k²) cost-matrix fill**: only unmatched splits (k ≪ n
    for similar trees) enter the cost-matrix fill and LAP.

Applied to `mutual_clustering_score`, `msd_score`, and `jaccard_score`.
Not applied to `shared_phylo_score` or `msi_score` (exact-match
detection is incorrect for SPI/MSI — see earlier bug fix).

**Files changed:** `src/pairwise_distances.cpp` (new
`find_exact_matches()` helper + refactored score functions).

**A/B benchmark (release build, same-process comparison):**

| Scenario               | ref (ms) | dev (ms) | Change   |
|------------------------|----------|----------|----------|
| CID 100×50-tip         | 70.1     | 16.0     | **−77%** |
| CID 40×200-tip         | 178.2    | 17.6     | **−90%** |
| MSD 100×50-tip         | 64.8     | 13.9     | **−79%** |
| MSD 40×200-tip         | 215.8    | 15.9     | **−93%** |
| PID 100×50-tip         | 122.9    | 95.1     | −23%     |
| IRF 100×50-tip         | 37.4     | 16.9     | −55%     |
| CID cross 20×30 50-tip | 24.8     | 3.9      | **−84%** |
| MSD cross 20×30 50-tip | 13.8     | 3.1      | **−77%** |
| LAPJV 400 (canary)     | 2.12     | 1.27     | neutral  |
| LAPJV 1999 (canary)    | 95.2     | 95.5     | neutral  |

Random trees (no exact matches): no regression (full LAP path
unchanged). Numerical accuracy: max \|ref − dev\| ≤ 5.7e-12.

Key insight: for MCMC posteriors and bootstrap replicates (the common
real-world case), most splits are shared between trees. The sort+merge
pre-scan reduces the effective LAP dimension from ~n to ~3, yielding
order-of-magnitude speedups.

### VTune hotspot profile — PID-focused (vtune-pid/, 2026-03-17)

**Workload**: PID only — similar 50/200-tip + random 50/200-tip, ~20 s
each (80 s total). Build flags: `-O2 -fopenmp -fno-omit-frame-pointer`
(pre-popcount spi_overlap).

| Function                        | CPU Time | %     | Notes                                                   |
|---------------------------------|----------|-------|---------------------------------------------------------|
| `lap`                           | 36.6s    | 45.8% | LAP solver — dominates PID since full LAP always solved |
| `spi_overlap`                   | 23.0s    | 28.7% | 4-pass boolean scan (pre-popcount)                      |
| `func@0x340aa1a00` (unresolved) | 14.9s    | 18.7% | Likely inlined `shared_phylo_score` cost-matrix fill    |
| `free_base`                     | 0.33s    | 0.4%  | Heap deallocation                                       |
| `malloc_base`                   | 0.28s    | 0.3%  | Heap allocation                                         |

**Key observations:** 1. **LAP at 45.8%** — with the SPI/MSI correctness
fix (always solving the full LAP, no exact-match shortcut), the LAP
solver is now the overwhelmingly dominant cost. The LAP has been
extensively optimized and offers no easy further gains. 2.
**`spi_overlap` at 28.7%** — a major target. Called n² times per pair to
fill the full cost matrix. The 4-pass boolean-scan approach was replaced
by popcount (below). 3. **malloc/free \< 1%** — CostMatrix pooling has
effectively eliminated allocation overhead. 4. The random 200-tip
workload dominates total time (~21s of 80s elapsed) because the full
197×197 LAP is solved for every pair with no shortcuts.

### Popcount-based `spi_overlap` (DONE, kept)

Replaced the 4-pass boolean-scan `spi_overlap` with a single-pass
popcount approach.

**Old approach** (4 sequential scans with early exit): 1. Scan bins for
`a & b` → if none found, A and B disjoint 2. Scan bins for `~a & b` → if
none found, B ⊆ A 3. Scan bins for `a & ~b` → if none found, A ⊆ B 4.
Scan bins for `~(a | b)` with last-bin mask → if none found, A ∪ B = all
tips 5. If all 4 regions populated → return 0

Each pass resets pointers and re-iterates, requiring 4× memory loads for
the common “return 0” case. For n_bins=4 (200-tip trees), each pass had
loop setup, pointer reset, and branch overhead.

**New approach** (single popcount pass):

``` cpp
int16 n_ab = 0;
for (int16 bin = 0; bin < n_bins; ++bin) {
  n_ab += TreeTools::count_bits(a_state[bin] & b_state[bin]);
}
// Derive all 4 Venn-diagram regions from n_ab, in_a, in_b, n_tips:
//   n_a_only  = in_a - n_ab
//   n_b_only  = in_b - n_ab
//   n_neither = n_tips - in_a - in_b + n_ab
// Case selection via 3-4 integer comparisons (no loops, no pointer resets)
```

Benefits: - **One pass** over `a_state[]` and `b_state[]` (each bin
loaded exactly once) - **No branches in scan loop** — just AND +
hardware POPCNT + ADD - **No pointer resets** between passes - **Handles
all n_bins uniformly** — same code for 1-bin and 4-bin cases -
`#include <TreeTools/SplitList.h>` added to `tree_distances.h` for
`TreeTools::count_bits` access (needed by `day_1985.cpp` which doesn’t
include SplitList.h directly)

**Files changed:** `src/tree_distances.h` (spi_overlap rewritten;
SplitList.h included).

**A/B benchmark (release build, same-process comparison via
`compare-ab.R`, both builds `-O2 -fopenmp`, no
`-fno-omit-frame-pointer`):**

| Scenario                | ref (ms) | dev (ms) | Change           |
|-------------------------|----------|----------|------------------|
| PID 100×50-tip          | 168      | 138      | **−18%**         |
| PID 40×200-tip          | 558      | 423      | **−24%**         |
| MSID 100×50-tip         | 164      | 155      | −5%              |
| MSID 40×200-tip         | 741      | 613      | **−17%**         |
| CID 100×50-tip (canary) | 28.6     | 30.1     | +5% (noise)      |
| CID 40×200-tip (canary) | 30.7     | 36.1     | +18% (alignment) |
| MSD 100×50-tip (canary) | 26.9     | 20.1     | −25% (alignment) |
| MSD 40×200-tip (canary) | 29.3     | 29.0     | neutral          |
| IRF (canary)            | —        | —        | neutral          |
| LAPJV (canary)          | —        | —        | neutral          |

Numerically exact (max \|ref − dev\| = 0 for all metrics).

The 200-tip PID improvement (−24%) is the strongest because n_bins=4
makes the single-pass popcount most advantageous vs the old 4-pass
approach.

**Canary alignment noise:** CID and MSD canaries show ±5–25% swings
between runs despite not calling `spi_overlap`. Three runs with
different flag combinations produced MSD 50-tip results of +14%, +3.5%,
and −25%. This is the same instruction-alignment sensitivity documented
for the LAP header split: changing the inline function body in
`tree_distances.h` shifts code layout in `pairwise_distances.o`,
affecting branch prediction and instruction fetch for neighbouring
functions. The swings cancel across metrics and scenarios, confirming
they are not systematic.

------------------------------------------------------------------------

## Combined A/B Benchmark: main vs dev (all optimizations)

**Baseline (ref):** main branch (OpenMP PR \#176 only). **Dev:** current
development branch (all optimizations including sort+merge +
MatchScratch pooling). Release build, same-process comparison via
`compare-ab.R`. Trees: `as.phylo(0:N, tipLabels = ...)` — similar trees
(high split overlap).

| Metric          | Scenario   | ref (ms) | dev (ms) | Change   |
|-----------------|------------|----------|----------|----------|
| CID             | 100×50-tip | 83.5     | 16.0     | **−81%** |
| CID             | 40×200-tip | 210.8    | 18.2     | **−91%** |
| MSD             | 100×50-tip | 82.4     | 12.3     | **−85%** |
| MSD             | 40×200-tip | 253.7    | 23.5     | **−91%** |
| PID             | 100×50-tip | 143      | 112      | −22%     |
| PID             | 40×200-tip | 408      | 358      | −12%     |
| MSID            | 100×50-tip | 157      | 114      | −27%     |
| MSID            | 40×200-tip | 474      | 462      | −3%      |
| IRF             | 100×50-tip | 52.1     | 14.3     | **−73%** |
| IRF             | 40×200-tip | 84.4     | 18.1     | **−79%** |
| CID cross 20×30 | 50-tip     | 20.6     | 3.9      | **−81%** |
| MSD cross 20×30 | 50-tip     | 17.5     | 3.8      | **−78%** |
| LAPJV 400       | (canary)   | 1.47     | 1.47     | neutral  |
| LAPJV 1999      | (canary)   | 120      | 112      | neutral  |

Correctness: all metrics max \|ref − dev\| ≤ 5.5e-12 (floating-point
level). LAPJV canary: neutral (no LAP regression).

CID, MSD, IRF benefit most (73–91%) from sort+merge exact-match
detection. PID and MSID show 3–27% from R-level fast paths, C++ batch
entropy, and CostMatrix pooling (exact-match detection is incorrect for
SPI/MSI — see earlier bug fix).

### ManyMany fast paths for \*Distance functions (DONE, kept)

Created `.FastManyManyPath()` helper in `tree_distance.R` that reuses
pre-computed Splits for both pairwise distance computation AND per-tree
entropy/info calculations. This complements `.FastDistPath()`
(all-pairs) for the cross-pairs case.

**Pattern**: When `tree1` and `tree2` are both collections with uniform
(identical) tip sets and `reportMatching=FALSE`: - Single `as.Splits()`
call per collection with unified tip labels - Pass both Splits to C++
batch cross-pairs function (`cpp_*_cross_pairs`) - Simultaneously
compute per-tree entropies via C++ batch entropy function - Compute
normalization info via `outer(info1, info2, "+")` - Return to calling
distance function (e.g., `ClusteringInfoDistance`)

**Applied to:** - `ClusteringInfoDistance` (tree_distance_info.R) -
`DifferentPhylogeneticInfo` (tree_distance_info.R) -
`MatchingSplitInfoDistance` (tree_distance_msi.R) - `InfoRobinsonFoulds`
(tree_distance_rf.R)

**Correctness verified**: max \|fast − slow\| ≤ 2.8e-14 (floating-point
level).

------------------------------------------------------------------------

## LinkingTo Header Exposure (`expose-lapjv` branch)

Extracted LAP and MCI C++ APIs into `inst/include/TreeDist/` headers so
that downstream packages (e.g., TreeSearch) can use
`LinkingTo: TreeDist`:

| Header                     | Content                                            |
|----------------------------|----------------------------------------------------|
| `types.h`                  | `cost`, `lap_dim`, `lap_row`, `lap_col`, constants |
| `cost_matrix.h`            | `CostMatrix` class (Rcpp-free)                     |
| `lap_scratch.h`            | `LapScratch` struct                                |
| `lap.h`                    | `lap()` declarations                               |
| `lap_impl.h`               | `lap()` implementation (include in exactly one TU) |
| `mutual_clustering.h`      | MCI declarations                                   |
| `mutual_clustering_impl.h` | MCI implementation (include in exactly one TU)     |

`src/lap.h` is now a thin wrapper that includes `<TreeDist/lap.h>` and
re-exports types to global scope.

### LAPJV codegen regression (diagnosed, characterised, accepted)

Including `lap_impl.h` in `lap.cpp` changed the TU context enough for
GCC 14’s register allocator to produce ~8% more instructions in the
Dijkstra hot loop, causing a consistent 20–25% regression on standalone
LAPJV (n ≥ 400).

**Root cause:** the installed-header version of `CostMatrix` (in
`inst/include/TreeDist/cost_matrix.h`) has a different method set than
main’s monolithic `src/lap.h` (extra methods like `setWithTranspose()`,
`dim8()`; missing test variants like `findRowSubminNaive`). This changes
GCC’s optimization heuristics for the entire TU, even though `lap()`
never calls the differing methods.

**Fix:** define `lap()` directly in `lap.cpp` (not via
`#include <TreeDist/lap_impl.h>`) with
`__attribute__((optimize("align-functions=64", "align-loops=16")))`. The
`lapjv()` wrapper fills the transposed buffer first (matching R’s
column-major storage) then untransposes — restoring the cache-friendly
pattern.

**VTune profiling (vtune-lap/, 2026-04-01):** Profiled with
`-O2 -g -fno-omit-frame-pointer` on the current branch (LAPJV 1999×1999
×30s + 400×400 ×30s + CID/MSD random 200-tip ×30s each).

| Function                              | CPU Time | %     | Notes                                               |
|---------------------------------------|----------|-------|-----------------------------------------------------|
| `TreeDist::CostMatrix::findRowSubmin` | 23.8s    | 19.8% | Inlined into lap() at -O2 (separate symbol with -g) |
| `TreeDist::lap`                       | 23.4s    | 19.5% | Dijkstra phase dominates                            |
| `msd_score`                           | 9.4s     | 7.8%  | Contains inlined lap()                              |
| `mutual_clustering_score`             | 9.3s     | 7.7%  | Contains inlined lap()                              |
| `lapjv` (R wrapper)                   | 6.0s     | 5.0%  | Matrix conversion overhead                          |

Source-line breakdown of `lap()`:

| Phase                       | Lines   | CPU Time (s) | % of lap() |
|-----------------------------|---------|--------------|------------|
| Dijkstra inner update loop  | 308–325 | ~16.5        | 72%        |
| Dijkstra min-scan           | 274–287 | ~3.5         | 15%        |
| Matrix fill (lapjv wrapper) | 61–65   | ~2.7         | 12%        |
| findRowSubmin (reduction)   | 173–184 | ~2.0         | 9%         |

**Optimisation attempts:**

1.  **`__restrict__` on Dijkstra pointers** (DONE, kept): Added
    restrict-qualified local pointers (`d_ptr`, `pred_ptr`, `cl_ptr`)
    for the augment-solution phase to eliminate potential alias reloads.
    Applied to both `lap.cpp` and `lap_impl.h`. Benchmark: no measurable
    change in regression magnitude, but correct optimisation.

2.  **`#pragma GCC optimize("O3")`** (attempted, reverted): Forced O3
    for the entire `lap.cpp` TU. Benchmark: no measurable improvement;
    same alignment pattern.

3.  **Per-object Makevars flags** (attempted, failed): R’s build system
    does not support per-target variable overrides in `Makevars.win`.

**Residual regression (confirmed, alignment-dependent):**

| LAPJV size | dev/ref ratio | Direction                             |
|------------|---------------|---------------------------------------|
| 400×400    | 0.81          | **Dev 19% faster** (alignment win)    |
| 1999×1999  | 1.08–1.13     | **Dev 8–13% slower** (alignment loss) |

The regression is alignment-sensitive and manifests differently at
different problem sizes. The same codegen change that makes 400×400
faster makes 1999×1999 slower. This is the classic instruction-alignment
lottery: the Dijkstra inner loop’s branch prediction and instruction
fetch are affected by code placement that varies with TU context.

**Conclusion:** The regression cannot be reliably eliminated without
matching main’s exact TU context (adding dead code back to the installed
header, which is unacceptable for CRAN). **Tree distance metrics are
completely unaffected** — `lap()` is called from
`pairwise_distances.cpp` (different TU context), and benchmarks confirm
neutral performance across all metrics and tree sizes.

**Maintenance note:** if the `lap()` algorithm changes, update BOTH
`src/lap.cpp` and `inst/include/TreeDist/lap_impl.h`. The duplication is
intentional — it preserves the TU context that was profiled and tuned.

------------------------------------------------------------------------

## Remaining Optimization Opportunities

- Sort+merge pre-scan for `rf_info_score`: **DONE** — replaced O(n²)
  loop with `find_exact_matches()` + O(k) info summation over matches.
  No LAP involved, so the entire per-pair cost drops to O(n log n) for
  the sort+merge step. All tests pass, numerically exact (max \|fast −
  slow\| = 0).
- `MatchScratch` allocation pooling: **DONE** — `find_exact_matches()`
  internal buffers (canonical forms, sort indices, match arrays) are now
  pooled via a `MatchScratch` struct, one per OpenMP thread. Eliminates
  per-pair heap allocation for sort+merge matching in all batch
  functions (CID, MSD, IRF, Jaccard all-pairs and cross-pairs).
  Release-build benchmark needed to measure impact.
- KendallColijn distance: **DONE** — replaced pure-R
  [`KCVector()`](https://ms609.github.io/TreeDist/dev/reference/KendallColijn.md)
  with C++ `cpp_kc_vector()` in `src/path_vector.cpp` (modelled on
  existing `path_vector()`). Per-tree vector speedup: **218× for 50-tip,
  226× for 200-tip** (debug build). Also optimised
  [`KendallColijn()`](https://ms609.github.io/TreeDist/dev/reference/KendallColijn.md)
  all-pairs and cross-pairs Euclidean distance using existing C++
  `pair_diff_euclidean()` and `vec_diff_euclidean()`. All existing tests
  pass.
- LAP inner loop: AVX2 SIMD intrinsics assessed and deemed
  **unpromising** — 64-bit cost values limit AVX2 to 4-wide (same as
  existing manual 4× unroll); indirect `col_list` indexing prevents
  efficient vectorisation; platform independence required for CRAN.
  Closed.
- SPR distance: algorithm-limited (see profiling section below). Closed.

------------------------------------------------------------------------

## SPR Distance Profiling (VTune, 2026-03-13)

### Workload

**vtune-spr-driver-v2.R**: 100 iterations each of: - Similar 150 ×
50-tip (11,175 pairs/iter) - Similar 60 × 200-tip (1,770 pairs/iter) -
Random 150 × 50-tip (11,175 pairs/iter)

**Total elapsed**: 342 s CPU / 339 s effective  
**TreeDist.dll contribution**: ~11 s (3.2% of total, rest is R overhead)

### Top TreeDist hotspots (by CPU time)

| Function                     | CPU Time | % of TreeDist.dll | Notes                                      |
|------------------------------|----------|-------------------|--------------------------------------------|
| `func@0x340ae4910` (unnamed) | 3.87s    | 36%               | Likely inlined merge of multiple functions |
| `reduce_trees`               | 1.35s    | 12.5%             | Tree reduction logic (core SPR algorithm)  |
| `keep_tip` (TreeTools)       | 0.36s    | 3.3%              | Tip filtering in tree manipulation         |
| `keep_and_reroot`            | 0.28s    | 2.6%              | Root manipulation during reduction         |
| `SplitList::SplitList`       | 0.28s    | 2.6%              | Rcpp conversion overhead                   |
| `calc_mismatch_size`         | 0.12s    | 1.1%              | Split mismatch computation (VTune target)  |

### Key findings

1.  **Tree reduction is the dominant bottleneck** (≈1.35s, 12.5%), not
    the split comparison (`calc_mismatch_size` is only 1.1%). This is
    architecturally fundamental to the deOliveira heuristic: the
    algorithm iteratively reduces trees until convergence.

2.  **The main unnamed function (3.87s, 36%)** likely represents inlined
    code from multiple smaller functions. The DLL was built without
    `-mpopcnt` in this profile, but that accounts for \< 1%.

3.  **The SPR algorithm** (`SPRDist` R → `keep_and_reduce` C++ →
    `reduce_trees` → `keep_and_reroot` → `TreeTools::keep_tip`)
    involves:

    - Finding the smallest mismatched split via `mismatch_size`
    - Computing disagreement split (XOR operation, R-level)
    - Selecting tips to keep based on disagreement
    - Recursively reducing both trees via `keep_and_reduce` →
      `ReduceTrees` (R) → `keep_and_reroot` → `reduce_trees` (C++)
    - Each reduction is a O(n) tree restructuring pass

4.  **The algorithm is inherently iterative**: SPRDist returns the
    number of reduction iterations performed plus 1. For the benchmark
    trees (similar 50-tip with high split overlap), the typical
    iteration count is **low** (~2–3 moves), suggesting the reduction
    loop terminates quickly — but the per-iteration cost is high.

### No obvious low-hanging fruit

- **Early termination heuristic**: Possible if we could detect “final
  state” earlier, but this requires algorithmic insight into tree
  topology constraints.
- **Parallel reduction**: The outer loop over iterations is inherently
  sequential (each iteration depends on the previous), so parallelism is
  not applicable.
- **SIMD in `reduce_trees`**: The algorithm is not data-parallel (tree
  pointers, linked-list manipulation, highly branching control flow).
  SIMD unlikely to help.
- **Vector vs pool**: The array allocations in `reduce_trees` are
  per-pair, with sizes fixed at call time (n_vert). Could pool them, but
  the complexity is high for likely \< 5% gain (similar to prior
  CostMatrix pooling).

### Conclusion

The SPR distance computation is fundamentally limited by the iterative
tree reduction algorithm. The per-pair cost scales with tree size (O(n))
and iteration count, which varies by tree similarity. Without changing
the algorithm (e.g., using an exact solver like TBRDist), significant
speedups are unlikely.

**Future optimization paths** (beyond scope of current dev cycle): 1.
**Investigate the Rogue method** (experimental, under development) — may
be faster 2. **Exact solver integration** via TBRDist R package (if
performance is critical for user) 3. **Batch SPR via MCMC** — amortize
tree reduction across proposals 4. **Algorithm literature** — post-2008
SPR heuristics may exist with better constants

**Recommendation**: Close SPR optimization; move to other metrics or
accept it as algorithmically constrained.
