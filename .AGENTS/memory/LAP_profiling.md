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

### Test suite fix

Added `tests/testthat/setup.R` that opens a null PDF device for the duration of
all test runs.  This suppresses bare `plot()` / `TreeDistPlot()` calls in tests
from appearing in the interactive graphics device, and prevents vdiffr snapshot
rendering from leaking to screen.  vdiffr opens its own `svglite` device on top
of the null device, so snapshot tests are unaffected.
