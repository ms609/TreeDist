# Is a C++ port of `KMeansPP` worthwhile?

Assessment after the O(n²) → O(n·dim) seeding refactor. Short answer: **the big
win is already banked; a C++ port buys ~1.3–1.5× at most and only at large n, so
it is not justified now.** Porting `stats::kmeans` is specifically *not*
worthwhile unless paired with an accelerated algorithm.

## Where the time goes (evidence)

Decomposition of `KMeansPP.matrix(x, k = 10, nstart = 10)`, dim = 36:

| n      | total | seeding loop (distance rows) | kmeans (×10) |
|--------|------:|-----------------------------:|-------------:|
| 6 400  | 1.85 s | 0.47 s (0.45 s)             | 1.26 s (~68%) |
| 20 000 | 5.31 s | 1.41 s (1.27 s)             | ~5.3 s (dominant) |
| 50 000 | 8.32 s | 3.64 s (3.25 s)             | 2.86 s (~34%) |

`Rprof` self-time at n = 50000 (the cleanest apportionment; the table above is
noisy because Hartigan–Wong's iteration count is data-dependent):

| symbol         | self % | what it is                          |
|----------------|-------:|-------------------------------------|
| `.Fortran`     | 49 %   | `kmeans` Hartigan–Wong core (compiled) |
| `.rowSums`     | 40 %   | distance-row computation (compiled C)  |
| `sample.int`   | 4 %    | D²-weighted draw (compiled)            |
| everything else| ~7 %   | kmeans internals (`aperm`/`sweep`/`colMeans`), `pmin.int`, R glue |

**~89 %+ of runtime is already in compiled kernels.** The R-level glue we wrote
(`.DistanceRow` closure, the `for` loops, `pmin.int`) is <1 % self-time. There is
no slow interpreted loop left to eliminate — the earlier refactor already did
that by deleting the O(n²) matrix build.

## What C++ could and could not buy

**1. Fused distance-row kernel (the only low-risk option).**
The current row costs ~3 transient n×dim allocations (`rep`, the subtraction, the
square) plus a separate `sqrt` pass. A fused C++ kernel would stream `x` once,
accumulate `Σ(x[i,j]−c[j])²` per row, and `sqrt` in place — no temporaries — and
could use the OpenMP already in `src/`. Realistic gain on the `.rowSums` portion
~2–3× single-thread (allocation + cache), more with threads.
*Amdahl ceiling:* distance rows are ~40 % at n=50000 but only ~24 % at n=6400, so
overall this is **~1.3–1.5× at large n, ~1.2× at n=6400**. Effort: ~50–100 lines
mirroring existing `src/` patterns. Reproducibility: same minor FP caveat as now
(fusing changes summation order slightly; would need the same identity re-check).

**2. Porting `stats::kmeans` to C++ — not justified.**
It is 49 % (n=50000) to ~68 % (n=6400) of the time, so it looks tempting, but it
is **already compiled, well-tuned Fortran**. A naïve C++ Lloyd/Hartigan–Wong
reimplementation would be a wash at best and likely slower (Hartigan–Wong
converges in fewer passes than textbook Lloyd). A genuine speedup requires an
*accelerated* algorithm — Elkan/Hamerly triangle-inequality bounds — which is
research-grade effort, must handle empty-cluster and convergence edge cases, and
**breaks exact reproducibility** (a C++ RNG cannot match R's `sample.int`, and
bounds-accelerated k-means returns equivalent-not-identical results). Payoff
maybe ~1.5–2× overall; cost and risk are out of proportion for a utility helper.

**3. A pure-R BLAS shortcut (mentioned for completeness, rejected).**
`dist²(i) = ‖xᵢ‖² + ‖c‖² − 2 xᵢ·c` turns the row into a BLAS `gemv` (fast, often
multithreaded) and avoids the `rep`. But it suffers catastrophic cancellation
for near-coincident points (needs `pmax(0, ·)`) and **changes the FP results**,
which would forfeit the bit-identity we just verified and could flip near-tie
draws. Deliberately not used.

## Recommendation

- **Do nothing in C++ for now.** The asymptotic + memory win (O(n²) → O(n·dim))
  is the headline prize and it is captured. The residual cost is ~90 % in
  already-compiled kernels, so there is no cheap interpreter overhead to remove.
- **If, and only if, a real workload profiles distance-row seeding as a
  bottleneck at very large n**, write the fused OpenMP distance-row kernel
  (option 1): best effort-to-reward ratio, reproducibility-preservable.
- **Do not port `stats::kmeans`.** If kmeans itself becomes the bottleneck at
  scale, the right lever is algorithmic (the already-cited scalable k-means‖
  seeding, Bahmani 2012, and/or Elkan/Hamerly refinement), not a transliteration
  of the existing algorithm into C++.

_Profiling driver inline in this directory's session notes; reproduce with
`Rprof` around `KMeansPP(matrix(rnorm(50000*36), ncol=36), k=10, nstart=10)`._
