# Pol & Escapa leaf instability as the `SPRDist(method = "Rogue")` selector

Research prototype for ms609/TreeDist#174 at 0dd26a81. No package code is changed. The scripts and results are in `dev/rogue_polescapa/`.

## Recommendation

**Adopt it as the selector, and keep MCI only to break ties** (variant `both-tied` below). Score each leaf by the Fitch length it adds, over both trees. Drop the highest-scoring leaf. When several leaves share the top score, score just those with the existing bounded MCI.

- **Accuracy.** It is better than the current selector on random pairs of 25 or more leaves, and much better at 50–150 leaves. Across the 60 extra large pairs it wins 52 and loses 1. Elsewhere it is equal, except on random pairs of 8–14 leaves, where it is slightly worse (below).
- **Cost.** The Fitch pass costs 25–125 µs per iteration in C++, against 1.5–99 ms for `spr_rogue_scores()`. Pure Fitch (`both-k0`) measured roughly 2× faster than `SPRDist` at HEAD at 30 leaves, ~10× at 80 and ~40× at 150. `both-tied` was not timed end to end, because the prototype calls the full `spr_rogue_scores()` on every tie. Its projected extra cost is scoring 2–3 tied candidates in 30–46% of iterations, about 0.1–6% of the current selector's MCI work.

Using it only as a pre-filter (top-k then MCI) is dominated. Small k gives no better accuracy than `both-tied` and costs k LAPs per iteration. Larger k drifts back towards MCI's worse choices at large n.

## Method

Each split of tree 2 is a binary character, scored on tree 1 with Fitch, and vice versa. A leaf's score is the number of steps saved when it is coded `?` in every character. Coding a leaf `?` is the same as deleting it, so one downpass and one uppass per character give every leaf's score: L − (cost of the tree minus that leaf's subtree). Each (character, leaf) term is 0 or 1.

- `selectors.R`: a harness (`RogueWith()`) that copies the depth-1 loop of `.SPRRogue()` with a pluggable leaf choice, the incremental Fitch (`FitchAmbiguous()`), and the variants.
- `validate.R`: checks two things. The incremental Fitch equals `TreeSearch::TreeLength()` with the leaf recoded `?` (40/40 pairs). The harness with the MCI scorer reproduces `SPRDist(method = "Rogue")` (60/60).
- `fitch_bitpar.cpp`: a bit-parallel C++ version (64 characters per word, two bit planes). It matches the R prototype on 40/40 pairs.
- `make_reference.R`, `evaluate.R`, `summarize.R`, `fitch_bench.R`: build the reference set, run the selectors, print `tables.md`, and time the C++ version.

Variants: `both` / `12` / `21` = both directions, tree 2's splits on tree 1, or the reverse. `k0` = pure Fitch, ties broken by leaf order (`-last`: by reverse order). `kK` = MCI among the K top-scoring leaves. `tied` = MCI among all leaves tied at the top score. `hyb20` = full MCI while the reduced trees have ≤ 20 leaves, Fitch above that. `mci-then-fitch` = current selector with Fitch breaking MCI ties. All runs are at depth 1; depth 2 is out of scope.

## Accuracy

Cells give the number of pairs over-estimated by 0 / +1 / +2 / ≥ +3 moves, with the mean in brackets. W/L counts pairs where the variant gives a lower / higher estimate than `mci`. Full tables are in `dev/rogue_polescapa/tables.md`.

**Against exact uSPR** (`TBRDist::USPRDist`). This covers 528 pairs:

- Random pairs of 8–14 and 15–20 leaves.
- Near pairs of 20–50 leaves, 2–6 or 7–10 leaf moves apart.
- The two fixtures from git history. Their trees regenerate exactly: 36/36 spot checks agree with the stored matrices.

25 of the 80 random 15–20-leaf pairs exceeded the 60 s limit for `USPRDist` and are excluded. That biases this set towards easier pairs.

No estimate fell below the exact distance. All 528 × 13 estimates are ≥ exact uSPR, as the leaf-dropping construction guarantees, and `summarize.R` asserts this.

| Selector | rand8-14 (120) | rand15-20 (55) | near20-50 (80) | near20-50x (59) | fix25-12spr (78) | fix32-16spr (136) | W/L |
|---|---|---|---|---|---|---|---|
| mci (current) | 111/9/0/0 (0.07) | 33/17/5/0 (0.49) | 80/0/0/0 (0) | 57/2/0/0 (0.03) | 58/20/0/0 (0.26) | 119/16/1/0 (0.13) | – |
| both-k0 | 96/23/1/0 (0.21) | 30/21/4/0 (0.53) | 80/0/0/0 (0) | 56/2/1/0 (0.07) | **78/0/0/0 (0)** | 116/20/0/0 (0.15) | 36/36 |
| 12-k0 | 89/28/3/0 (0.28) | 23/23/9/0 (0.75) | 79/1/0/0 (0.01) | 50/5/3/1 (0.24) | 78/0/0/0 (0) | 114/20/2/0 (0.18) | 29/54 |
| both-k3 | 108/11/1/0 (0.11) | 34/18/3/0 (0.44) | 80/0/0/0 (0) | 58/0/1/0 (0.03) | 65/13/0/0 (0.17) | 118/18/0/0 (0.13) | 22/14 |
| both-k5 | 111/9/0/0 (0.07) | 33/18/4/0 (0.47) | 80/0/0/0 (0) | 59/0/0/0 (0) | 65/10/3/0 (0.21) | 119/17/0/0 (0.12) | 19/11 |
| **both-tied** | 99/21/0/0 (0.17) | 32/20/3/0 (0.47) | 80/0/0/0 (0) | 58/0/1/0 (0.03) | **78/0/0/0 (0)** | 118/18/0/0 (0.13) | 34/24 |
| hyb20-k5 | 111/9/0/0 (0.07) | 33/17/5/0 (0.49) | 80/0/0/0 (0) | 57/2/0/0 (0.03) | 58/20/0/0 (0.26) | 118/17/1/0 (0.14) | 0/1 |

**Against best-known** (the minimum over all variants), for 60 extra random pairs (`evaluate.R extra`), where exact uSPR is out of reach:

| Selector | rand50 (30) | rand80 (20) | rand150 (10) | W/L |
|---|---|---|---|---|
| mci (current) | 2/8/11/9 (1.97) | 0/4/2/14 (3.35) | 0/0/2/8 (4.60) | – |
| both-k0 | 14/9/6/1 (0.80) | 8/6/5/1 (0.95) | 4/4/1/1 (0.90) | 50/5 |
| both-k5 | 13/11/5/1 (0.80) | 9/5/3/3 (1.05) | 3/4/0/3 (1.40) | 45/5 |
| **both-tied** | **20/7/2/1 (0.47)** | **13/3/3/1 (0.60)** | 5/3/1/1 (0.90) | **52/1** |
| hyb20-k0 | 17/6/6/1 (0.70) | 9/4/7/0 (0.90) | 2/3/4/1 (1.40) | 49/4 |

The `dev/rogue_bound` workloads tell the same story. On rand30, rand80 and rand150, `both-tied` has a mean excess of 0.42 / 1.00 / 2.50, against 0.90 / 2.67 / 4.75 for `mci`, and wins 30 / loses 6. On near80, near150, mid80 and mid150, every variant that uses both directions has zero excess, as `mci` does.

Observations:
- **Both directions are needed.** One direction alone is worse than `mci` almost everywhere. The two directions see different conflict.
- **Leaf order does not drive the result.** Breaking pure-Fitch ties by reverse order (`both-k0-last`) moves the means only within noise: 0.21 → 0.23 and 0.42 → 0.57.
- **MCI helps only as a tie-breaker.** Widening the MCI choice (k = 3 → 5 → 10) moves the large-tree results back towards `mci`'s worse ones. Applying MCI to small reduced trees inside large runs (`hyb20`) also costs accuracy.
- **The one regression is small random trees** (8–14 leaves). `both-tied` over-estimates by 1 on 21 of 120 pairs, against 9 of 120 for `mci`. At that size both selectors are cheap, so a production version could run full MCI on input pairs of ≤ ~14 leaves. That cut-off was not tested; `hyb20` shows the cut-off must apply to the input size, not to the reduced trees inside a run.

## Fitch tie sets

The top score is unique in most iterations. The share of iterations with a tie (tie set > 1) is:

| Workload | Share tied |
|---|---|
| Near pairs | 2–11% |
| Random pairs, 8–30 leaves | 27–33% |
| Random pairs, 50–80 leaves | 32–40% |
| Random pairs, 150 leaves | 43–46% |

A tied set averages 2.0–3.2 leaves and has ≥ 5 leaves in at most 8% of iterations. Over a whole pair, `both-tied` runs MCI on 0.1–6% of the candidates that `mci` scores (1.4% at 150 leaves).

## Cost

All timings are on this machine, with another session's benchmark running at the same time, so they are indicative only. "Per iteration" is the median for one scoring step on the first pair of each workload.

| Leaves | Fitch kernel, C++, both directions | + R glue (`as.Splits`, `RootTree`) | `spr_rogue_scores()` | `SPRDist` per pair at HEAD | Harness + C++ pure Fitch, per pair |
|---|---|---|---|---|---|
| 30 | 25–30 µs | 0.3 ms | 1.5 ms | 16 ms | 7.8 ms |
| 80 | 50–55 µs | 0.5 ms | 17 ms | 0.40 s | 34 ms |
| 150 | 125–170 µs | 0.8–1.0 ms | 99–108 ms | 4.0–4.4 s | 108 ms |

On `mid80` / `mid150`, the harness takes 5.8 / 18 ms per pair, against 51 / 510 ms for `SPRDist`.

The kernel times include packing the characters from an R logical matrix (n² scalar reads), which a native version would skip. The kernel is O(nodes × ⌈splits/64⌉) per direction: one pass up the tree and one down, with a few AND/OR/popcount operations per word. In TreeDist the characters are the other tree's split bitsets, transposed, so a C++ version needs no phyDat and no new data structures. The harness-plus-Fitch figures above are dominated by R glue and `spr_keep_and_reduce()`, so a native version would be faster still.

For `both-tied`, add the MCI tie-break on top. It scores 2–3 candidates in 30–46% of iterations, and those candidates could skip the parent-price LAP. That adds a few percent of the current MCI work.

## Caveats

- Exact references stop at 32 leaves. Larger sets compare against best-known, which is a relative measure.
- The 150-leaf conclusions rest on 14 pairs (10 extra plus 4 from the workloads). The direction of the effect is consistent: `mci` was never best.
- Depth-2 pair scoring and `sprProxy` alternatives are not covered.
- Tie-breaks differ from the current selector; the user confirmed bit-identical tie-breaks are not required.

## Next steps

1. **Implement `both-tied` in `src/spr.cpp`.** The Fitch kernel takes the two trees' split lists and returns per-leaf scores. `spr_rogue_scores()` would then score only the tied leaves. This touches the same file as the incremental-fill work, so it needs to follow that work.
2. **Test a full-MCI fallback** for input pairs of ≤ ~14 leaves, and decide whether the regression on small random trees matters.
3. **Re-run the depth-2 comparison** once the depth-1 selector has changed.
