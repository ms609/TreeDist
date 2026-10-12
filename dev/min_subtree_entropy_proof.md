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

## 1. The complete tree is optimal while `f` is non-decreasing

**Theorem 1.** Let `f` be non-decreasing and concave on `{2, ..., s}`.
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
`[0, n/2]` and stays at 1 above `n/2`.  Because `Hb` has zero slope at
`1/2`, `ĝ` is non-decreasing and concave on the reals, and hence on the
integers.

**Theorem 2.**  For `n >= 2` and `1 <= s <= n` (with `z = 0` when `s = n`),

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

**Corollary 2 (splits).**  Map an optimal `R` back to `T`.  If `H`'s first
merge is with an item of weight `X < s`, the root of `T` splits
`X | s - X`.  If `H` is merged last, with all `s` units, the root of `T`
splits as the last unit merge did.  So every Huffman tie-break gives an
optimal `k`.

*Conjecture (tie-completeness).*  The converse also holds: every optimal `k`
arises from some Huffman tie-break.  It is checked exactly against the DP
(optimal `k` taken within 1e-11) for every `n <= 4000` and every `s`,
and at sampled `n` up to 131 074.  Every non-candidate `k` was at least
6.9e-7 worse for `n <= 4000`, so there are no hidden exact ties there.  It
matters only when distinct `k` tie in exact arithmetic.  If it failed, the
table could differ from an exhaustive scan by rounding of a tied value, and
`cost` would still be exact to rounding.

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
3 candidates in every tested `n`, found in
`O(log s)` steps.  `huffman_splits()` implements this.

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

For `s > n/2`, bit-identity with the exhaustive scan is not proven.  It
needs the tie-completeness conjecture and a floating-point gap.  The
smallest non-candidate gap is about `2 / (n² ln 2)`.  The crude error bound
`n · 7e-15` from Section 1 overtakes that gap near `n ≈ 7 · 10^4`, while
the measured error (about 1e-13) would not until about `5 · 10^6`.  So, for
`s > n/2`, bit-identity rests on the conjecture plus the experiments below.
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

## 6. Prior art

- Glassey & Karp (1976, *SIAM J. Appl. Math.* 31:368–378, Thm 1 and §4)
  show that Huffman trees minimize the sum of any non-decreasing concave
  function of internal-node weights.  With equal weights, the Huffman tree
  is the complete tree.  Theorem 2 rests on this result.
- Cleary, Fischer & St. John (2025, *Bull. Math. Biol.*,
  doi:10.1007/s11538-025-01522-1, Thm 4.1 and Lemma 4.17) prove the
  strictly increasing, strictly concave case of Theorem 1, with uniqueness.

The proof of Theorem 1 above, via `min(c, t)`, is an independent route.
Theorem 2 appears to be new: it reduces a symmetric, non-monotone split cost
to a monotone one by unrooting at the heavy leaf and capping at `n/2`.  The
same argument applies to any concave symmetric `f`.  The literature search
was limited.  Hu, Kleitman & Tamaki (1979), Parker (1980) and Golumbic
(1976) have not been read directly.
