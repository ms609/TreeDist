# Optimal subtree shapes for `.MinSubtreeEntropy()`

`min_subtree_entropy_impl()` (`src/nni_max_step.cpp`) fills, for a tree of
`n` leaves,

    cost[1] = 0,   cost[s] = f(s) + min_{1 <= k <= s/2} (cost[k] + cost[s - k]),

with `f(c) = Hb(c / n)`, the binary entropy (bits) of a `c | n - c` split.
`cost[s]` is the least value of `F(T) = sum over internal nodes v of f(c_v)`
over rooted binary trees `T` with `s` leaves, where `c_v` is the number of
leaves below `v`.

Write `C_s` for the complete ("heap-shaped") binary tree with `s` leaves:
a perfect tree on `2^L <= s` leaves whose `s - 2^L` leftmost leaves are each
replaced by a cherry.  Its root splits `k* = max(P/2, s - P) | s - k*`, where
`P` is the largest power of two below `s` (`complete_split()`).

**Attribution.**  None of the structural results below is new; see Section 6.
- Theorem 1 is Glassey & Karp (1976, §4).  Cleary, Fischer & St. John (2025,
  Thm 4.1) give the strictly concave case with uniqueness.
- Theorem 2 is a special case of Goubko (2016, Thm 3 with Lemmas 19–20).
- Theorem 3 follows from Goubko's equality clause (Thm 3, Lemma 8).
- Our own contributions are:
  - the heavy-leaf reduction of the recurrence above to Goubko's setting;
  - the short capping proof of Theorem 2 via Glassey & Karp, and the
    independent flip proof of Theorem 3;
  - the `O(log s)` enumeration of every optimal split (Section 3), hence the
    `O(n log n)` table.

## 1. The complete tree is optimal while `f` is non-decreasing

**Theorem 1** (Glassey & Karp 1976, §4)**.** Let `f` be non-decreasing and concave on `{2, ..., s}`.
Then `F(C_s) <= F(T)` for every rooted binary tree `T` with `s` leaves.

*Reduction.*  With `D_u = f(u+1) - f(u)` and `d_u = D_u - D_{u+1} >= 0`
(`d_{s-1} = D_{s-1} >= 0`),

    f(c) = const + sum_{u=2}^{s-1} d_u * min(c, u + 1),   2 <= c <= s.

Every tree has `s - 1` internal nodes, so it suffices to show, for every
integer `t >= 2`, that `C_s` minimizes

    Phi_t(T) = sum_v min(c_v, t).

*Tools.*  `E(d)`, the least external path length of a tree with `d` leaves,
is attained by `C_d` (Knuth) and is convex and piecewise linear:
`E(x) = x(j + 2) - 2^(j+1)` on `[2^j, 2^(j+1)]`.  For any tree,
`sum_v c_v` equals its external path length.

*Lower bound.*  If `s <= t`, `Phi_t(T) = sum_v c_v >= E(s) = Phi_t(C_s)`.
Otherwise let `U = {v : c_v > t}`; it contains the root and is closed under
taking ancestors.  The maximal clades with `c <= t` ("pendants") number
`m = |U| + 1`, have sizes `d_i` summing to `s`, and each contributes its
own external path length `>= E(d_i)`.  So, by Jensen with the convex
interpolant `Ê` of `E`,

    Phi_t(T) >= t(m - 1) + sum_i E(d_i) >= psi(m) := t(m - 1) + m Ê(s / m).

`psi` is convex in real `m` (a perspective).  Where `s/m` lies in
`[2^j, 2^(j+1)]`, `psi'(m) = t - 2^(j+1)`.  Let `2^J <= t < 2^(J+1)`.  The
slope is negative for `j >= J` and non-negative for `j < J`, so the integer
minimum is at `m = floor(s / 2^J)` or `ceil(s / 2^J)`.  Writing
`s = a 2^J + b`, `0 <= b < 2^J`, gives

    Phi_t(T) >= Lambda_t(s) := t(a - 1) + aJ 2^J + b(J + 1) + min(b, t - 2^J).

*Attainment.*  In `C_s` (with `s = 2^L + r`) the `N = 2^(L-J)` nodes at
height `J` have sizes `2^J + clamp(r - p 2^J, 0, 2^J)`, `p = 0..N-1`: `q`
nodes of `2^(J+1)`, at most one of `2^J + b`, and the rest `2^J`, where
`a = N + q`.  The `N - 1` nodes above them exceed `t` and score `t` each.  A
node of `2^(J+1) > t` scores `t + 2E(2^J)`.  A node of `2^J` scores `E(2^J)`.
The partial node scores `E(2^J + b)` if `2^J + b <= t`.  Otherwise it scores
`t + (J - 1)2^J + b(J + 1)`, because its daughters lie in `[2^(J-1), 2^J]`.
Summing gives exactly `Lambda_t(s)`.  ∎

(The identity `Phi_t(C_s) = Lambda_t(s)` was also checked exhaustively for
all `s, t <= 700`.)

**Corollary 1 (uniqueness and gap).**  If `f` is strictly concave, then
`F(T) - F(C_s) >= min_{2 <= u <= s-2} d_u` unless `T` has the same multiset
of clade sizes as `C_s`.  The counts `#{v : c_v >= t}` determine the
multiset.  `Phi_2 = 2(s - 1)` and `Phi_s - Phi_(s-1) = 1` hold for every
tree.  So if the multisets differ, some `Phi_t` with `3 <= t <= s - 1`
differs by at least 1.  The restriction matters: at odd `n` with
`s = ceil(n/2)`, `d_(s-1) = 0`.  In the DP, the candidate `k` stands for the tree
`(C_k, C_(s-k))`.  Its largest non-root clade is `s - k`, so its multiset
equals that of `C_s` only if `k = k*`.

**Application.**  `Hb(c/n)` is non-decreasing on `c <= ceil(n/2)`, and
strictly concave with `|f''(c)| = 1 / (ln 2 · c (n - c)) >= 4 / (n² ln 2)`.
Hence, for `s <= ceil(n/2)`, `k*` is the unique minimizer.  Every other `k`
is worse by at least `4 / (n² ln 2)` (5.8e-8 at n = 10^4).

Floating-point error in `cost[]` is at most about `n · 7e-15`:
- `|error in h| <= ~6u log2 n`, given a `log2` correct to 1 ulp, summed
  over at most `n/2` nodes;
- plus addition rounding `<= 2u s log2 s`.

That is below half the gap for `n` up to about 45 000.  So, for
`s <= ceil(n/2)`, `split[s] = k*` and
`cost[s] = h[s] + (cost[k*] + cost[s - k*])` reproduce the exhaustive scan
bit for bit.  The implementation uses exactly these operations.

Above `n ≈ 45 000` the guarantee is empirical only.  The gap bound is loose:
the measured smallest gap is about `0.43 / n`, and the measured worst error
is 3.4e-11 at `n = 45 000`.  Bit-identity was also checked directly
(Section 4).

## 2. Every size: Huffman trees with one heavy leaf

Let `z = n - s` and `ĝ(w) = Hb(min(w, n/2) / n)`.  `ĝ` equals `f` on
`[0, n/2]` and stays at 1 above `n/2`.  `ĝ` is non-decreasing and concave
on the reals, and hence on the integers: on `[0, n/2]` its slopes are
non-negative and non-increasing, and beyond that its slope is 0.  The zero
slope of `Hb` at `1/2` is not needed; see *Generalization* below.

**Theorem 2** (a special case of Goubko 2016, Thm 3)**.**  For `n >= 2` and `1 <= s <= n` (with `z = 0` when `s = n`),

    cost[s] = Huf(s) - 2 + ĝ(z),

where `Huf(s)` is the least value of `sum over internal v of ĝ(w_v)` over
rooted binary trees `R` with leaf weights `{z, 1, ..., 1}` (`s` ones).  Here
`w_v` is the total leaf weight below `v`.  So, by Glassey & Karp (Section 6),
`Huf(s)` is attained by Huffman's algorithm.

*Proof.*

- **Unrooting.** Join a leaf `H` of weight `z` by an edge to the root of a
  tree `T` counted by `cost[s]`.  This gives an unrooted binary tree `U`,
  and every `U` with leaf `H` arises from exactly one `T`.  Write `g~(e)` for
  `Hb(w / n)` at an edge splitting the weights `w | n - w`; `g~` is
  symmetric.  Then `cost(T) = sum over e in U of g~(e) - s ĝ(1)`, because
  the root of `T` corresponds to `H`'s edge, with split `z | s`.
- **Rooting.** Root `U` by subdividing an edge `e0`, giving `R`.  Each
  non-root node `u` has `ĝ(w_u) >= g~(edge above u)`, since `ĝ >= g~`.  The
  two root children together score `ĝ(w1) + ĝ(w2) = g~(e0) + 1`, because the
  heavier has `w >= n/2`.  So
  `sum over non-root u of ĝ(w_u) >= sum over e of g~(e) + 1`.
- **Equality.** Pick `e0` incident to a centroid: a vertex all of whose
  branches weigh `<= n/2`.  The centroid may be the leaf `H` itself, when
  `z > n/2` or `s = 1`; then `e0` is `H`'s edge.  Every non-root node except
  the centroid side's root child weighs `<= n/2`, so the bound holds with
  equality.  Every rooted `R` arises this way from some `U`.
- **Counting.** The sum over non-root nodes equals
  `sum over internal v of ĝ(w_v) - ĝ(n) + s ĝ(1) + ĝ(z)`, and `ĝ(n) = 1`.
  Minimizing over `U` gives `Huf(s) = cost[s] + 2 - ĝ(z)`.  ∎

All Huffman tie-breaks give the same multiset of merged weights: the two
smallest values at each step are fixed, whichever equal items are taken.
So every Huffman tree is optimal.

Theorem 1 is the special case `z >= s`.  Then `H` is merged last, and the
units form a Huffman tree on equal weights, which is the complete tree.

*Generalization.*  Nothing above uses the shape of `Hb` beyond monotonicity
and concavity on `[0, n/2]`.  Let the split cost be `h(min(c, n - c))`, with
`h` non-decreasing and concave on `[0, n/2]`.  Equivalently, the cost is any
symmetric concave function of `c` on `[0, n]`.  Put
`ĥ(w) = h(min(w, n/2))`, which is concave for the same reason as `ĝ`.  The
same proof then gives

    cost[s] = Huf_ĥ(s) - 2 h(n/2) + ĥ(z).

Theorem 3 below additionally needs `h` strictly concave on `[0, n/2]`.  It
fails for linear `h`: at `n = 8`, `s = 6` the optimal `k` are `{1, 2, 3}`,
but Huffman gives only `{2}`.  So, for every strictly concave `h`, the
optimal-split set is the same as for `Hb`.  A third review confirmed both
statements in exact arithmetic for `h` = sqrt, log1p and `w^0.7`, all with
non-zero slope at `n/2`, for `n <= 40`.

**Corollary 2 (splits).**  Map an optimal `R` back to `T`.  If `H`'s first
merge is with an item of weight `X < s`, the root of `T` splits
`X | s - X`.  If `H` is merged last, with all `s` units, the root of `T`
splits as the last unit merge did.  So every Huffman tie-break gives an
optimal `k`.

**Theorem 3 (tie-completeness).**  Conversely, for `2 <= s <= n`, every
optimal `k` arises from some Huffman tree, i.e. some tie-break of Huffman's
algorithm.

*Proof via Goubko (2016).*  Root each unrooted tree at a centroid vertex,
or at the leaf's neighbour if the centroid is a leaf (Goubko's Lemma 19).
These are Goubko's proper trees, and Huffman trees are among them
(Lemma 20).
- In a proper tree every non-root internal vertex has subordinate weight
  `<= n/2`, where `ĝ` is strictly increasing and strictly concave.
- Goubko's Thm 2 (compare Forst & Thorup's Thm 5.5) makes the Huffman weight
  vector weakly majorize every other proper tree's.  By strictness, the cost
  is then strictly lower unless the sorted vectors coincide.
- If they coincide, Goubko's Lemma 8 makes the tree a Huffman tree.

So every optimal unrooted tree is an unrooted Huffman tree, and `k` depends
only on the unrooted tree.

Two caveats apply.  Goubko's Lemma 11 asserts strictness for every
increasing concave `χ`, which is false as stated (take `χ` linear).  It holds
for strictly increasing, strictly concave `χ`, which covers `Hb(x/n)` on
`[0, n/2]`.  Goubko also requires positive pendant weights, so `s = n`
(`z = 0`) is covered only by the proof below and by the exact checks
(Section 4).  ∎

*Independent proof via Forst & Thorup (2005).*  The hypothesis `s > n/2` is
not needed; it was verified unnecessary for every `2 <= s <= n <= 300`
(Section 4).  For `s <= ceil(n/2)`, Corollary 1 also gives the result
directly.  Call `R` optimal if it attains `Huf(s)`.

- **Optimal `R` have a short heavy chain.**  In *Rooting*, a non-root node
  that is not a root child has `ĝ(w_u) > g~` exactly when `w_u > n/2`.  So
  if such a node exists, `R` scores strictly more than its own `U`'s bound,
  which is at least the optimum.  In an optimal `R`, only the root and one
  root child can weigh more than `n/2`.
- **`k` depends only on `U`.**  `k` is the smaller of the two non-`H`
  branches at `H`'s neighbour in `U`.  Re-rooting `R` without changing `U`
  leaves `k` alone.
- **Flips.**  Forst & Thorup (2005, Thm 5.3) transform any weighted tree
  into a Huffman tree by a finite sequence of strict H-flips.  The theorem
  as stated allows any strict H-flip.  The restriction to Lemma 5.2 flips
  comes from its proof: Lemma 5.2 supplies a flip of that form whenever the
  tree is not Huffman, and strict flips strictly decrease the weight vector
  in Schur order (Thm 5.5) on a finite set.  So we may assume every flip has
  that form.  Their Thm 5.5 and Cor 5.7 show that
  `G_ĝ = sum ĝ(internal weights)` never increases under an H-flip, since
  `ĝ` is non-decreasing and concave.  Start from an optimal `R`.  The end
  is Huffman, hence optimal, so every flip has zero cost and every
  intermediate tree is optimal.
- **Zero-cost Lemma 5.2 flips only re-root.**  In their notation an H-flip
  swaps `U` and `V`, at distance `l` from `V` to the common ancestor `V_l
  = U_{l+p}`.  The proof of Thm 5.5 changes the internal weights as:
  pairs `{u_i, v_i} -> {ũ_i, ṽ_i}` (`i < l`), with equal sums and `ũ_i`
  the smallest of the four; and singles `u_i -> ũ_i <= u_i` (`i >= l`).
  Each pair term is `<= 0` by concavity and each single term `<= 0` by
  monotonicity.  A zero-cost flip therefore has every term zero.
  Lemma 5.2 gives either `l = 1` (the uncle case) or (5.2.1).  (5.2.1) gives
  `v < u` and `u' < v'`, and also `u' <= u` and `v <= v'`.  It does not give
  `u < v'`; the argument below uses only `v < u` and `u' < v'`.  In the second case with `l >= 2`, pair
  `i = 1` goes from `{u + u', v + v'}` to `{v + u', u + v'}`.  `v + u'` is
  strictly the smallest of these four, so the pair changes.  (A general
  strict H-flip need not change a pair, so the restriction to Lemma 5.2
  matters.  For example, take `n = 7`, `s = 6`, and let `H` have weight 1.
  `T = (((H,(1,2)),(3,4)),(5,6))` has `k = 2`.  Flipping `U = (1,2)` with
  `V` = leaf 5 gives `(((H,5),(3,4)),((1,2),6))`, with `k = 1`.  Internal
  weights go from `{7,5,3,2,2,2}` to `{7,4,3,2,2,2}`, and both 5 and 4
  exceed `n/2`, so the cost is unchanged.  Pair `i = 1` goes `{3,2} -> {2,3}`,
  because `u' = v' = 1`, which violates `u' < v'`.)  A changed pair keeps `G_ĝ` only if `ĝ` is affine on
  `[ũ_i, ṽ_i]`.  That needs `ũ_i >= n/2`, as `Hb` is strictly concave on
  `[0, 1/2]`.  But `u_i + v_i <= n`, since `U_i` and `V_i` are disjoint, so
  `ũ_i = u_i = v_i = n/2`, which contradicts the change.
  So `l = 1`, `V` is the brother of `U_p`, and `p >= 1`.  The weights
  `u_1, ..., u_p` each fall by `u - v > 0`.  Zero cost then needs
  `u_1 - (u - v) >= n/2`.  So `U_1, ..., U_p` weigh more than `n/2`, are
  non-root, and lie below `U_{p+1}`.  By the first point, `p = 1`, `U_1` is
  the heavy root child and `V` the other root child.  In `U`, node `U_1`
  touches `U`, `U'` and `V`.  The flip moves the root from edge `U_1–V` to
  edge `U_1–U` and leaves `U` unchanged.
- So the optimal `R` and the Huffman tree at the end of the sequence share
  `U`, and hence `k`.  ∎

Section 3 lists `k` for every Huffman tree.  Ties among units leave `H`'s
partner weight unchanged.  Ties involving `H` are the `b = z` cases it
enumerates.  So the candidates are exactly the
optimal `k` in exact arithmetic.  For `n <= 4000`, every non-candidate `k`
was at least 6.9e-7 worse.

## 3. Huffman on `{z, 1, ..., 1}` in `O(log s)`

Unit-only Huffman runs level by level.  At each level it holds `count`
items of weight `reg = 2^j` and at most one other item of weight
`odd in (reg/2, 2 reg)`.  Merges come in this order:

1. `(odd, reg)` if `odd < reg`;
2. the `(reg, reg)` pairs;
3. if one `reg` is left, `(reg, odd)`, or `(reg, 2 reg)` if there is no
   `odd`.

The second-smallest weight `b` never decreases.  `H` first merges at the
first unit merge `(a, b)` with `b > z`, and takes `a` as its partner.  A
merge with `b = z` is a tie that `H` may win, so its `a` is also a
candidate.  If no merge has `b >= z`, `H` joins last.  That gives at most
3 candidates in every tested `n` (`n <= 2000` exhaustively), found in
`O(log s)` steps.  The `O(n log n)` bound needs only `O(log s)` candidates,
at most one per level step.  The `<= 3` bound is empirical.
`huffman_splits()` implements this.

Huffman coding by runs of equal weights, in time that depends on the number
of distinct weights rather than on the number of leaves, is due to Moffat &
Turpin (1998).  Here there are two distinct weights, 1 and `z`.
`huffman_splits()` specializes that idea to track `H`'s partner.  (Moffat &
Turpin was not read; it is cited on the editor review's advice.)

This also explains the shapes seen earlier:

- the spine of power-of-two pendants is `H` absorbing `reg`-sized items;
- the complete tree persists while `H` joins last;
- the near-root runs of constant `k`, or constant `s - k`, come from
  constant partner weights.

## 4. Implementation and verification

`min_subtree_entropy_impl()` uses `complete_split(s)` for `s <= ceil(n/2)`
(Corollary 1).  For larger `s` it evaluates the Huffman candidates in
increasing `k` and keeps the first floating-point minimum, as the scan
would.  The table build is `O(n log n)`.

For `s > n/2`, Theorem 3 makes the candidates exactly the optimal `k`.
Bit-identity with the exhaustive scan also needs a floating-point gap,
which is not proven.  The
smallest non-candidate gap is about `2 / (n² ln 2)`.  The crude error bound
`n · 7e-15` from Section 1 overtakes that gap near `n ≈ 7 · 10^4`, while
the measured error (about 1e-13) would not until about `5 · 10^6`.  So, for
`s > n/2`, bit-identity rests on Theorem 3 plus the experiments below.
`cost` is exact to rounding either way.  The branch-and-bound scan, its SIMD
helpers and the earlier pruning lemma are no longer needed.

`cost` and `split` are identical (`!=` on doubles) to a plain exhaustive DP
for:
- every `n <= 3000`;
- every 173rd `n` up to 30 000;
- `n = 100 003`.

An independent review found the same, comparing every `n <= 4000` and
sampled `n` up to 131 074 against a plain exhaustive DP.  The review also
checked Theorem 2's formula against an exhaustive search over trees with a
heavy leaf, for every `n <= 450`, without relying on Glassey & Karp.  It
compared `huffman_splits()` with a heap-based Huffman under all tie-breaks
over 1.17 million `(s, z)` cases.  No differences were found.

A second review checked Theorem 3 in high-precision (`Decimal`) arithmetic.
`Decimal` is not exact for logarithms.  For
`s <= 8`, it enumerated every rooted tree on `{z, 1^s}`:
- every optimal tree reached a Huffman tree by zero-cost re-rootings;
- the optimal `k` matched `huffman_splits()` and every Huffman tie-break.

The optimal `k` sets also matched for `n <= 160`.

A third review (2026-10) used exact arithmetic.  It stored each
`n ln2 Hb(c/n) = n ln n - c ln c - (n - c) ln(n - c)` as an integer vector of
prime exponents.  Logs of primes are linearly independent over the
rationals, so two values tie exactly when their vectors are equal.  Strict
inequalities were decided in floating point, falling back to 60-digit
`mpmath` below a gap of 1e-6; no near-tie below 1e-45 occurred.  For every
`2 <= n <= 300` and `1 <= s <= n` it confirmed:
- Theorem 2, against an exhaustive DP over all rooted trees on `{z, 1^s}`;
- Corollary 2 and Theorem 3, including `s <= n/2`;
- that the candidates from `huffman_splits()` are exactly the brute-force
  tie-break set.

It also found the level description identical to heap-based Huffman for
`s < 3000`, and at most `1.9 log2 s` level steps.

They are identical to the previous build for every `n <= 6000` and at
`n = 200 000`.  The block bounds passed on to `nni_max_step()` are
unchanged.

Table build times, C++ only:

| n | before | after |
|---|---|---|
| 300 | 17.5 µs | 5.2 µs |
| 1000 | 96 µs | 27 µs |
| 5000 | 1.14 ms | 0.35 ms |
| 20 000 | 6.9 ms | 1.6 ms |

## 5. History

An earlier version of this note (reviewed) proved only Theorem 1 and a
pruning lemma for `s >= 2n/3`.  Theorem 2 subsumes both.

Theorems 2 and 3 were at first believed new.  A wider search in 2026-10
found both in Goubko (2016); Section 6 now records this.  The results stand
unchanged, but the attribution was corrected.  The stand-alone write-up was
dropped in favour of a lemma in a paper on `NNIMaxStep()`.

## 6. Prior art

Read directly:
- Glassey & Karp (1976), Forst & Thorup (2005), Parker (1980) and Golumbic
  (1976), the first two again in 2026-10 against the cited result numbers;
- Goubko (2016), from arXiv:1502.01216 v2.  Theorem and lemma numbers below
  follow arXiv v2; the published MATCH version was not seen.

**Summary.**
| Here | Source | Status |
|---|---|---|
| Thm 1 | Glassey & Karp 1976, §4 (Cleary et al. 2025, Thm 4.1, strict case) | known |
| Thm 2 | Goubko 2016, Thm 3 + Lemmas 19–20 | known; our capping proof is shorter |
| Thm 3 | Goubko 2016, Thm 3 equality clause + Lemma 8 (strict `χ`) | known; our flip proof is independent |
| Reduction of `cost[]` to one heavy leaf | — | not found in prior art |
| Section 3 enumeration, `O(n log n)` table | (Moffat & Turpin 1998 for run-length Huffman) | not found in prior art; elementary |

**Goubko (2016)**, "Minimizing Wiener index for vertex-weighted trees with
given weight and degree sequences", *MATCH Commun. Math. Comput. Chem.*
75(1):3–27; arXiv:1502.01216.  This is the closest prior work.
- *Setting.*  Trees with prescribed vertex weights `µ(v) >= 0` and degrees
  `d(v)`.  Pendant weights must be positive, and weights must be
  degree-monotone among internal vertices (Def. 4).  Our case has internal
  vertices of degree 3 and weight 0, and leaves `{z, 1^s}`.  This fits
  whenever `z > 0`, i.e. `s < n`.
- *Generalized Huffman* (Def. 6).  Repeatedly join the lightest items into
  a star whose centre is the lightest free internal vertex.  For degree-3
  internals this is binary Huffman, with the last three items joined at one
  root vertex: the unrooted form of a binary Huffman tree.  Tie-breaks give
  the collection `T_H`.
- *Lemma 19.*  Every tree has a vertex whose deletion leaves components of
  weight `<= µ̄/2` (a weighted centroid).  Rooting there, or at a leaf
  centroid's neighbour, gives a "proper" directed tree.
- *Lemma 20.*  Every directed Huffman tree is proper.
- *Thm 3.*  For `χ` concave and increasing on `[0, µ̄/2]`, a Huffman tree
  minimizes `sum over non-root v of χ(f(v))` among proper trees, where
  `f(v)` is the weight subordinate to `v`.  Equality holds only for Huffman
  trees.  Every non-root internal vertex of a proper tree has
  `f(v) <= µ̄/2`.  So, up to a constant from the pendant edges, this sum
  equals `sum over edges e of χ(min(w_e, µ̄ - w_e))`.  With
  `χ = Hb(·/n)` this is Theorem 2 here, once the heavy leaf is attached.
- *Thm 1* (`χ(x) = x(µ̄ - x)`, the vertex-weighted Wiener index of Klavžar &
  Gutman 1997): the minimizers are exactly `T_H`.
- *Lemma 8.*  A tree with the same subordinate-weight vector as a Huffman
  tree is a Huffman tree.
- *Caveat.*  Lemma 11, quoted from majorization texts, claims strictness for
  every increasing concave `χ`.  That is false for linear `χ`.  It holds for
  strictly increasing, strictly concave `χ`, which is all Theorem 3 here
  needs.
- His Thm 2, that the Huffman subordinate-weight vector weakly majorizes
  every other tree's, follows Zhang et al. (2008, unit weights).  It is the
  analogue of Forst & Thorup's Thm 5.5.

**Glassey & Karp (1976)**, "On the optimality of Huffman trees", *SIAM J.
Appl. Math.* 31:368–378, doi:10.1137/0131030:
- Their Theorem 1 says a tree minimizes `sum f(W_j)` over internal weights
  for *every* non-decreasing concave `f` if and only if it is a Huffman tree
  (H-tree).  This holds for any non-negative leaf weights and any tie-break,
  including our weight-0 leaf.  They need `f` on `[0, ∞)`; extend `f`
  linearly between integers and as a constant beyond `s`.
- They characterize H-trees as the trees with a lexicographically minimal
  sorted weight sequence.
- Section 4 solves the recurrence of Theorem 1 above,
  `F(n) = f(n) + min_j [F(j) + F(n - j)]`, `1 <= j <= n - 1`, with `f`
  concave and non-decreasing.  The minimizing `j` is the unique power of two
  with `n/3 <= j < 2n/3`.  With `k <= s/2` as in our recurrence, one part of
  the minimizing split is that power of two.  For `s = 7` the power of two
  is 4, so `k = 3`.  This matches `complete_split()` (checked for all
  `s < 2 · 10^5`).
- Their Lemma 2 proves the majorization step with the functions
  `min(x, s_q)`, the same reduction as Section 1.

So Theorem 1 is theirs, and Section 1 is only a self-contained re-proof.

Hong's converse, their Theorem 2, assumes `f ∈ C²[0, ∞)`.  For any such `f`
that is not non-decreasing and concave, some weights make every H-tree
suboptimal.  `Hb(c/n)` is not non-decreasing; its non-monotone case is the
one that matters, so H-trees cannot be applied to it directly.  Capping at
`n/2` after unrooting at the heavy leaf (Theorem 2 here) gets around this;
Goubko gets around it with proper rootings.

**Forst & Thorup (2005)**, "What costs are minimized by Huffman trees?",
*SIAM J. Discrete Math.* 19:46–68, doi:10.1137/S0895480103421725:
- They characterize Huffman trees by subtree-weight inequalities (Thm 2.3).
- Every tree reaches a Huffman tree by strict H-flips (Thm 5.3).
- The internal-weight multiset is a strict Huffman cost under Schur order
  (Thm 5.5).
- Hence a strictly increasing, strictly concave `g` is minimized only at
  Huffman trees (Cor 5.7; the converse is Note 5.6).  "Increasing" there is
  weak.

`ĝ` is flat above `n/2`, so Cor 5.7's uniqueness does not apply to rooted
trees.  The flip proof of Theorem 3 uses Thm 5.3 and the proof of Thm 5.5
instead.  The step showing that zero-cost Lemma 5.2 flips only re-root is
ours.  Its conclusion also follows from Goubko's equality clause.

**Parker (1980)**, "Conditions for optimality of the Huffman algorithm",
*SIAM J. Comput.* 9:470–489, doi:10.1137/0209035, extends Glassey & Karp in
two ways: to "quasilinear" weight-combination functions, and to Schur-concave
costs of the internal weights (Theorems 5 and 9).  It does not treat
symmetric or non-monotone node costs.  An erratum exists, not read: *SIAM J.
Comput.* 27:317 (1998), doi:10.1137/S0097539797328550.

**Cleary, Fischer & St. John (2025)**, "The GFB tree and tree imbalance
indices", *Bull. Math. Biol.* 87(10):145, doi:10.1007/s11538-025-01522-1.
Their Thm 4.1 rediscovers the strict case of Theorem 1 with uniqueness,
without citing Glassey & Karp.

**Golumbic (1976)**, *IEEE Trans. Comput.* C-25:1164–1167: merging with
cost `1 + max` of the children.  Not relevant.

**Related, lower relevance or not read:**
- *Cai, Cao, Li & Wang (2018)*, "On distances in vertex-weighted trees",
  *Appl. Math. Comput.* 333:435–442, doi:10.1016/j.amc.2018.03.117.  **Not
  read** (closed access).  The abstract claims "a simpler approach (than
  that of Goubko's) to obtain a stronger result" for a given degree
  sequence.  It may state the general-`χ` or the tie case, so check it
  before citing either as Goubko's alone.
- *Klavžar & Gutman (1997)*, *Discrete Appl. Math.* 80:73–81,
  doi:10.1016/S0166-218X(97)00070-X: defines the vertex-weighted Wiener
  index.
- *Zhang, Xiang, Xu & Pan (2008)*, "The Wiener index of trees with given
  degree sequences", *MATCH* 60:623–644, and *Wang (2008)*, *Discrete Appl.
  Math.* 156:2647–2654, doi:10.1016/j.dam.2007.11.005: greedy trees minimize
  the Wiener index for a given degree sequence.  This is the unit-weight
  case.  Zhang et al. was not read.
- *Schmuck, Wagner & Wang (2012)*, *MATCH* 68:273–292: greedy trees are
  extremal for many Wiener-type (distance-based) invariants.
- *Moffat & Turpin (1998)*, "Efficient construction of minimum-redundancy
  codes for large alphabets", *IEEE Trans. Inf. Theory* 44:1650–1657,
  doi:10.1109/18.681345: run-length Huffman construction.  Not read.
- *Hu, Kleitman & Tamaki (1979)*, "Binary trees optimum under various
  criteria", *SIAM J. Appl. Math.* 37:246–256, doi:10.1137/0137015.  Not
  obtained.  The abstract and citing papers concern alphabetic
  (order-preserving) and minimax trees; provisionally low relevance.
- Ruled out after checking:
  - Knuth (1982, algebraic Huffman);
  - Fujiwara & Jacobs (2014, leaf-depth cost functions);
  - Gallager (1978, sibling property; linear cost);
  - Golomb (1980) and Schwartz (1964) on Huffman tie-breaks;
  - Dasgupta (2016) hierarchical-clustering cost: it depends on a similarity
    graph, and on unit cliques all trees tie;
  - Lemant et al. (2022), whose entropy balance index is not `Hb(c/n)`.

Full search records (queries, ledgers, referee reports) are in the
2026-10 review artifacts, not in this repository.
