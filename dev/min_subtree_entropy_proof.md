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
(Section 3).

## 2. Above `n/2`

`f` is symmetric, `f(c) = f(n - c)`, and decreases beyond `n/2`.  Clades of
more than `n/2` leaves form a chain (spine) from the root.  By Theorem 1,
every other clade can be taken to be complete.

**Lemma 2 (pruning).**  Let `n/2 < s`, `k <= s/2`, and `a = complete_split(k)`
(the smaller daughter of `C_k`, so `k - a` is the larger).  Suppose
`n - s + a <= n/2` (true for all `k` once `s >= 2n/3`).  If
`k - a > n - s`, then `k` is not optimal for `s`; `a` does strictly better.

*Proof.*  `cost[k] = f(k) + cost[a] + cost[k - a]` (Theorem 1, `k <= n/2`).
Re-hang the daughter `k - a` beside `s - k`.  The tree with root split
`a | s - a`, whose `s - a` side splits `k - a | s - k`, scores
`f(s - a) - f(k)` relative to the tree using `k`.  Here
`f(s - a) = f(n - s + a)`, and `n - s + a < k <= n/2`, where `f` is strictly
increasing.  So `cost[a] + cost[s - a]` is less than `cost[k] + cost[s - k]`,
by at least one step of `f` below `n/2` (about `2.9 / n²`).  ∎

The hypothesis `n - s + a <= n/2` is in fact implied, because
`n - s + a < k <= n/2`.  For `s < 2n/3`, however, `m > s/2`, so the cap
below has no effect there.

*Floating point.*  Both candidates reuse the same computed `cost[a]`,
`cost[k - a]` and `cost[s - k]`.  Rounding is monotone, and `cost[s - a]`
is at most the rounded `h[s - a] + (cost[k - a] + cost[s - k])`.  The
comparison therefore reduces to `h[n - s + a] < h[k]`, a true gap of at
least `2.9 / n²`, against about four roundings of size `<= u·n`.  That is
safe up to `n ≈ 1.9 · 10^5`.  Because `a < k`, even a tie resolves to the
smaller index, so the first argmin is never a pruned `k`.

`k - complete_split(k)` is non-decreasing in `k`, so the surviving `k` form
an initial segment.  With `m = n - s` and `Q` the largest power of two
`<= m`, the segment is `k <= m + Q` (or `k = 1` if `m = 0`).  Near the root
only a handful of candidates remain.

## 3. Implementation and verification

`min_subtree_entropy_impl()` uses Corollary 1 for `s <= ceil(n/2)`.  For
`s >= 2n/3` it caps `k` at `m + Q` (Lemma 2).  Otherwise it keeps the
existing branch-and-bound scan.  Against the previous build, `cost` and
`split` are identical (`!=` on doubles) for every `n <= 3000`, every 97th
`n` up to 20 000, and `n = 50 000`.  An independent review compared the new
build against a plain exhaustive DP and found it identical for:
- every `n <= 6000`;
- every 37th `n` up to 20 000;
- spot values up to `n = 131 073`.

The same review re-derived each step of Theorem 1 and Lemma 2, and checked
them numerically.  The table build is 1.6–1.9× faster
(n = 300: 25 → 16 µs; 1000: 121 → 69 µs; 5000: 2.6 → 1.4 ms).

The table is still O(n²), because `n/2 < s < 2n/3` and the lower part of
`s >= 2n/3` are scanned.

## 4. Open: an O(n) build

Observed in the data, but not proven:

- `C_s` stays optimal up to a threshold between `0.715 n` and `0.795 n`,
  which varies irregularly with `n` (n <= 1200, 2000, 3001, 5000).
- Above it, the optimal tree is a spine.  The pendants, from the root down,
  are one subtree of arbitrary size and then powers of two, increasing toward
  a centroid node.  One daughter of that centroid node has `P/2` leaves.
  Example, n = 2000, s = 1976: pendants 24, 32, 128, 256, 512, then a
  perfect 1024.
- The larger daughter `s - split[s]` is non-decreasing in `s`.
- Distinct splits can tie exactly (n = 19, s = 15: `k = 3` and `k = 4`).

Two proven necessary conditions also hold on the spine.  Pendants increase
away from the root, because swapping two neighbours changes one spine clade
in `f`'s increasing range.  Each pendant exceeds the spine weight two steps
nearer the root, because otherwise merging the two pendants costs no more.
So the spine has `O(log n)` pendants.  The open difficulty is that, unlike in
Section 1, no single tree is optimal for every `min(c, n - c, t)` at once.
So the shape above `n/2` depends on the particular values of `Hb`, and a
closed form will need a function-specific argument.

The complete-tree cost is not convex, so Monge/SMAWK speed-ups don't apply
either.

## 5. Prior art

Theorem 1 is known.

- Glassey & Karp (1976, *SIAM J. Appl. Math.* 31:368–378, Thm 1 and §4)
  show that Huffman trees minimize the sum of any non-decreasing concave
  function of internal-node weights.  With equal weights, the Huffman tree is
  the complete tree.
- Cleary, Fischer & St. John (2025, *Bull. Math. Biol.*,
  doi:10.1007/s11538-025-01522-1, Thm 4.1 and Lemma 4.17) prove the
  strictly increasing, strictly concave case, with uniqueness, for the
  complete ("GFB") tree.

The proof above, via `min(c, t)`, is an independent route to the same
result.  I found no treatment of the non-monotone, symmetric case
`f(c) = Hb(c/n)`, so Lemma 2 and the open problem in Section 4 appear to be
new.  That search was limited, and Hu, Kleitman & Tamaki (1979), Parker
(1980) and Golumbic (1976) were not read directly.
