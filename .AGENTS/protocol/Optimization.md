# Benchmarking

All benchmarks live in `benchmark/` and use the `bench` package.

```r
source("benchmark/_init.R")   # loads TreeTools, TreeDist; defines Benchmark()
```

`Benchmark(expr)` wraps `bench::mark()` with `min_iterations = 3` and `time_unit = "us"`.
When run non-interactively (e.g. CI) results are serialised to `.bench.Rds` files.
`_compare_results.R` compares PR results against `main` and fails the build on a
regression > 4 %.

## ⚠ Benchmarking protocol — always use a release install

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

### Preferred workflow: A/B comparison (`compare-ab.R`)

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

### Legacy workflow: `bench_release()` / `bench_compare()`

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

## Existing benchmark scripts

| Script | What it tests |
|--------|---------------|
| `bench-tree-distances.R` | CID, PID, RF, MCI on 100×50-tip and 40×200-tip tree sets |
| `bench-LAPJV.R` | LAPJV on 40×40, 400×400, 1999×1999 uniform random matrices |
| `bench-PathDist.R` | PathDist on 6×182-tip trees |
| `bench-MCI-openmp.R` | MCI/CID OpenMP scaling on 100×50, 40×200, and 200×200-tip tree sets |

To add a new benchmark, create `benchmark/bench-<topic>.R` following the existing pattern
and add it to `_run_benchmarks.R`.

---

# Profiling with VTune

TODO this file needs streamlining.  Use the information below to infer a suitable vtune protocol. Check it works. Then REPLACE the below with your improved, tested protocol.

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


### VTune hotspot collection workflow (updated)

R's default `DLLFLAGS` includes `-s` (strip-all), which removes all symbols from the
DLL and makes VTune show addresses like `func@0x340a9ca80` instead of function names.
To produce a profiling build with full symbols:

1. Ensure `src/Makevars.win` exists with `-fno-omit-frame-pointer` in `PKG_CXXFLAGS`.
2. Override `DLLFLAGS` via `MAKEFLAGS` when installing to `.vtune-lib/`:
   ```r
   old_mf <- Sys.getenv("MAKEFLAGS")
   Sys.setenv(MAKEFLAGS = "DLLFLAGS=-static-libgcc")
   install.packages(".", lib = ".vtune-lib", repos = NULL, type = "source",
                    INSTALL_opts = "--no-multiarch")
   Sys.setenv(MAKEFLAGS = old_mf)
   ```
3. Run the VTune hotspot collection:
   ```
   vtune -collect hotspots -result-dir vtune-current \
         -- Rscript benchmark/vtune-driver.R
   ```
   (`benchmark/vtune-driver.R` exercises CID and PID workloads for ~30 s each.)

**Do not** add `DLLFLAGS` to `src/Makevars.win` — it has no effect there because
R's `Makeconf` sets `DLLFLAGS` after the package Makevars is parsed.  The `MAKEFLAGS`
environment variable overrides it correctly at the `make` command level.

After profiling, remove `-fno-omit-frame-pointer` from `src/Makevars.win`.