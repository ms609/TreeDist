# `KMeansPP.matrix` scaling: O(n²) → O(n·dim)

Results for [`kmeanspp_matrix_scaling.R`](kmeanspp_matrix_scaling.R), which compares
the old `KMeansPP.matrix()` (materialized `as.matrix(dist(x))`, O(n²) time and
memory) against the new on-the-fly seeding (per-centre distance rows, O(n·dim)).

**Environment.** R-devel (2026-06-02 r90096), x86_64-w64-mingw32, Windows 10;
TreeDist installed from source to a per-agent library (compiled, not `load_all`).
Timing via `system.time()`; peak memory via `gc(reset = TRUE)` then `gc()`
"max used" (low hundreds of MB are dominated by session baseline + transient
high-water, so treat them as indicative, not exact).

## 1. Reproducibility — results are unchanged

Under a fixed seed, the new seeding selects the **same centres** as the old
materialized-`dist` implementation, so cluster assignments and `tot.withinss`
are `identical()`:

| dataset                | cluster identical | `tot.withinss` identical |
|------------------------|:-----------------:|:------------------------:|
| n=2000, dim=5,  k=10   | ✅                | ✅ (5006.88)             |
| n=3000, dim=36, k=20   | ✅                | ✅ (93622.6)             |

The k-means++ RNG draw sequence (number and order of `sample.int` calls) is
preserved; only the numeric value of `min_d` can differ, by ≤ a few ULPs, which
did not flip any prob-weighted draw on the data tested. This is also asserted in
the unit suite (`tests/testthat/test-kmeanspp.R`, "matches materialized distance
matrix"), whose oracle is the literal pre-refactor body.

## 2. Timing + peak memory (k = 10, nstart = 10, dim = 36)

| n      | old: time / peak mem      | new: time / peak mem  | speed-up | old dense matrix |
|--------|---------------------------|-----------------------|----------|------------------|
| 6 400  | 4.6–6.0 s / ~781 MB       | ~1.8 s / ~142 MB      | ~2.6–3.2× | 0.33 GB         |
| 20 000 | *(not run)*               | ~5.2 s / ~442 MB      | —         | ~3.2 GB         |
| 50 000 | *(not run)*               | ~8.1 s / ~248 MB      | —         | ~20 GB          |

- **n = 6 400** (the diagnosed real-data size): the new method is **~2.6–3.2×
  faster and uses ~5.5× less peak memory**. The old ~5 s is dominated by the
  O(n²) `as.matrix(dist())` build (~3.3 s); the new method skips it entirely.
- **n = 20 000 / 50 000**: the old method must allocate an n×n double matrix —
  ~3.2 GB and ~20 GB respectively (peak ~1.5× that during `as.matrix(dist())`) —
  so it is effectively unusable. The new method completes in seconds with peak
  memory in the low hundreds of MB. (The new peaks are not strictly monotonic in
  n because they are `gc()` high-water readings; the substantive point is they
  stay **sub-GB** at every size while the old requirement grows as n².)

The old method was **not run** at n ≥ 20 000: a ~20 GB allocation would thrash
the page file rather than fail cleanly, so its cost is reported analytically
(n² × 8 bytes).

## 3. Asymptotics

| method               | seeding time            | memory      |
|----------------------|-------------------------|-------------|
| old `KMeansPP.matrix`| O(nstart · (n² + k·n))  | **O(n²)**   |
| new `KMeansPP.matrix`| O(nstart · k · n · dim) | **O(n·dim)**|

The `kmeans()` refinement step (coordinate space) is unchanged and shared by
both. `KMeansPP.dist` remains O(n²) by construction — it clusters the n rows of
an n×n distance matrix as n-dimensional vectors — and is unchanged here beyond
reusing the already-materialized matrix in its `kmeans()` call (avoids one
redundant n×n coercion).

## Reproduce

```sh
R CMD INSTALL --library=.dev-kmpp .
Rscript -e ".libPaths(c('.dev-kmpp', .libPaths())); source('tests/benchmark/kmeanspp_matrix_scaling.R')"
```
