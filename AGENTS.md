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

### Existing benchmark scripts

| Script | What it tests |
|--------|---------------|
| `bench-tree-distances.R` | CID, PID, RF, MCI on 100×50-tip and 40×200-tip tree sets |
| `bench-LAPJV.R` | LAPJV on 40×40, 400×400, 1999×1999 uniform random matrices |
| `bench-PathDist.R` | PathDist on 6×182-tip trees |

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

## Known Optimization Opportunities / TODOs

- `information.h` line 19: comment suggests considering increasing `MAX_FACTORIAL_LOOKUP`
  beyond 8192 or falling back to runtime calculation for larger values.
- `information.h` lines 120–122: code duplication flagged for consolidation.
- LAP inner loop: the 4× manual unroll works well; investigate whether AVX2 SIMD
  intrinsics (`_mm256_*`) could replace it on modern hardware.
- `CostMatrix` transpose: currently maintained as a full second copy; a cache-oblivious
  blocking scheme (already partially implemented via `BLOCK_SIZE`) could reduce memory
  bandwidth further.
- SPR distance (`spr.cpp`, `spr_lookup.cpp`): the algorithm is relatively recent
  (v2.8.0); profiling under VTune may reveal further hot spots.
- Parallelism: `parallel.R` exposes a parallel interface at the R level; the C++ layer
  does not yet use OpenMP. Pairwise distance matrices are an embarrassingly parallel
  workload.
- Large-tree path (`int32` migration, v2.12.0 dev): ensure new code paths are as
  optimized as the original `int16` paths.
