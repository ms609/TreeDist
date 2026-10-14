# Tighter MCI upper bounds for `SPRDist(method = "Rogue")`

Research on TreeDist at e07f2e9 (head of #211, since merged). Scratch only, nothing pushed. The patches are in `patches/`, next to this file.

## Recommendation

Use **parent-dual prices**, both as a bound and as a warm start (`patches/parent_prices_warm_start.patch`, +550/−20 lines, applies cleanly to e07f2e9).

1. At the start of each `spr_rogue_scores()` call, solve the MCI assignment once for the unreduced pair. Keep its optimal column duals v, shifted to [0, BIG].
2. Give each candidate's reduced tree-2 split the dual of the parent split it restricts to. Where splits merge, take the least dual.
3. After the cost matrix is filled, and only if the row bound has not already exited, compute `LB = Σ_j v_j + Σ_i min_j (c_ij − v_j)` in one exact int64 pass. If the implied MCI bound is below `min_score`, return it.
4. Otherwise run JV warm-started from v (`lap_warm()`).

**Lemma.** For any square cost matrix c and any vector v, LB(v) ≤ OPT.

*Proof.* For any permutation σ, `Σ_i c_iσ(i) = Σ_i (c_iσ(i) − v_σ(i)) + Σ_j v_j ≥ Σ_i min_j (c_ij − v_j) + Σ_j v_j`. ∎

So the bound is valid whatever the prices are: poor prices only make it looser. Validity therefore does not depend on the leaf mapping or on the warm start.

## Results

**Share of currently wasted LAPs rejected** (LAPs that run only to end below `min_score`):

| Bound | rand30 | rand80 | rand150 | near80 | near150 | Median excess at 80 tips (bits) |
|---|---|---|---|---|---|---|
| Entropy / size profile | 0% | 0% | 0% | 0% | 0% | 18–20 |
| Column reduction | 25% | 12% | 12% | 60% | 66% | 0.39 |
| Cold JV duals after ARR (the #174 idea) | 63% | 48% | 39% | 78% | 75% | 0.18 |
| Previous candidate's duals | 58% | 64% | 67% | 67% | 71% | 0.11 |
| **Parent column duals** | **82%** | **87%** | **88%** | **73%** | **86%** | **0.045** |

For scale, the median gap `min_score − MCI` on these LAPs is 0.14–0.20 bits.

**Wall-clock** (seconds, minimum of interleaved runs; output checksums identical across builds):

| Workload | Baseline | Cold post-ARR exit | Parent prices + warm start |
|---|---|---|---|
| 30 tips, random | 0.423 | 0.380 (−10%) | 0.343 (−19%) |
| 30 tips, depth 2 | 3.659 | 3.293 (−10%) | 3.131 (−14%) |
| 80 tips, random | 3.997 | 3.442 (−14%) | 1.944 (−51%) |
| 150 tips, random | 18.209 | 15.795 (−13%) | 6.396 (−65%) |
| 80 tips, 15 moves | 0.242 | 0.222 (−8%) | 0.219 (−10%) |
| 150 tips, 30 moves | 1.307 | 1.300 (−1%) | 1.001 (−23%) |
| 80 tips, near-identical | 0.048 | 0.045 | 0.047 (noise) |
| 150 tips, near-identical | 0.053 | 0.049 | 0.055 (noise) |

The LAPs that remain run 2.0×, 3.4× and 4.1× faster at 30, 80 and 150 tips. The LAP's share of runtime falls from 52–66% to 4–5%. The new overheads are the parent LAP (2–3%), the price mapping (2–3%) and the LB pass (1–4%).

**Correctness:**
- SPRDist gives identical results on 345 pair runs: 210 pairs of 8–160 leaves at depth 1, and 135 at depth 2.
- All 650,676 candidate scores (from 5,391 calls) were re-scored by both builds. The minima, argmin sets and near-minimum scores were identical, and the largest difference was 7e-15.
- Warm and cold LAPs differ by at most 1 cost unit (3.5e-15 bits), which is within the solver's existing ARR tolerance.

## Ideas that did not work

- **Size-profile and within-tree compatibility bounds** are about 100× too loose. They cannot see the overlaps between the two trees, which are what the cells measure.
- **Tightening the bound as rows are filled** would skip at most 1.5% of the fill, because leaving even one or two rows unfilled adds more slack than the margin.
- **Coordinate ascent** from column reduction stalls after one round.
- **Extra refinement on top of parent prices** (another coordinate round, or a warm ARR check) rejects 5–10 more points, but each costs another O(d²) pass, so the net effect is within ±1%.

## The earlier #174 verdict

#174 concluded that a post-ARR dual exit gave no measurable gain. That was measured when MCI was under 0.2% of runtime, before scoring moved to C++. Re-tested now, a cold post-ARR exit is exact and saves 10–14% on random pairs. Parent prices beat it everywhere.

## Next steps

1. **Productionize the patch.**
   - Add tests for the priced exit, the warm path, the mapping-failure branch, unequal split counts and the d > SL_MAX_SPLITS/2 guard.
   - Add a NEWS.md entry and a version bump.
   - Compute parent prices lazily, so near-identical pairs do not pay for the parent LAP.
2. **Speed up the cost-matrix fill**, which now takes 37–55% of runtime.
   - 97.6% of candidates drop a single leaf. A child cell's quadrant counts are the parent cell's counts minus that leaf's quadrant, so each cell costs O(1).
   - This would also skip KeepAndReduce and GetSplits (about 30% of runtime) for those candidates.
   - It is the largest remaining lever.
