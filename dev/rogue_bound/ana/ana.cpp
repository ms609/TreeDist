// Offline evaluation of upper bounds on the bounded MCI calls dumped by the
// instrumented TreeDist build (spr.cpp / instr_dump()).
//
// Usage: ana <dump.bin> <out.csv>
#include <x86intrin.h>
static unsigned long long g_init[2], g_arr[2], g_aug[2], g_n[2];
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>



#include "lap_hooked.h"

using namespace TreeDist;
using i128 = __int128;

static double lg2tab[4096];
static void InitLg2() {
  lg2tab[0] = 0;
  for (int i = 1; i < 4096; ++i) lg2tab[i] = std::log2(double(i));
}

struct SplitSet {
  int n_tip = 0, n_split = 0, nw = 0;
  std::vector<uint64_t> bits;
  std::vector<int> in;
  const uint64_t *row(int i) const { return &bits[size_t(i) * nw]; }
};

static bool ReadSplits(FILE *f, int n_tip, SplitSet &s) {
  int hdr[2];
  if (fread(hdr, sizeof(int), 2, f) != 2) return false;
  s.n_tip = n_tip;
  s.n_split = hdr[0];
  const int nc = hdr[1];
  s.nw = (n_tip + 63) / 64;
  s.bits.assign(size_t(s.n_split) * s.nw, 0);
  s.in.assign(s.n_split, 0);
  std::vector<uint8_t> buf(nc);
  for (int r = 0; r != s.n_split; ++r) {
    if (fread(buf.data(), 1, nc, f) != size_t(nc)) return false;
    for (int c = 0; c != nc; ++c) {
      s.bits[size_t(r) * s.nw + c / 8] |= uint64_t(buf[c]) << (8 * (c % 8));
    }
    int cnt = 0;
    for (int w = 0; w != s.nw; ++w) cnt += __builtin_popcountll(s.bits[size_t(r) * s.nw + w]);
    s.in[r] = cnt;
  }
  return true;
}

// ---- Weights exactly as mutual_clustering<true>() ----
static inline double IcMatching(int a, int b, int n) {
  return (a + b) * lg2tab[n] - a * lg2tab[a] - b * lg2tab[b];
}

static inline cost CellCost(int na, int nb, int a_and_b, int n) {
  const int nA = n - na, nB = n - nb;
  const int a_and_B = na - a_and_b;
  const int A_and_b = nb - a_and_b;
  const int A_and_B = nA - A_and_b;
  if (a_and_b == A_and_b && a_and_b == a_and_B && a_and_b == A_and_B) return BIG;
  const double lg2_n = lg2tab[n];
  const double offset_a = lg2_n - lg2tab[na];
  const double offset_A = lg2_n - lg2tab[nA];
  const double lg2_nb = lg2tab[nb], lg2_nB = lg2tab[nB];
  const double ic_sum =
    a_and_b * (lg2tab[a_and_b] + offset_a - lg2_nb) +
    a_and_B * (lg2tab[a_and_B] + offset_a - lg2_nB) +
    A_and_b * (lg2tab[A_and_b] + offset_A - lg2_nb) +
    A_and_B * (lg2tab[A_and_B] + offset_A - lg2_nB);
  const double max_over_tips = double(BIG) / n;
  return BIG - cost(ic_sum * max_over_tips);
}

// Reduced square LAP problem after removing exact matches.
struct Problem {
  int n_tip = 0, d = 0;
  double exact_part = 0; // exact-match score, already divided by n
  std::vector<int> rows, cols; // original split indices; -1 = padding
  CostMatrix C;
  double entropy_bound = 0;
};

static void Build(const SplitSet &a, const SplitSet &b, Problem &p) {
  const int n = a.n_tip;
  p.n_tip = n;
  std::vector<int> a_match(a.n_split, -1), b_match(b.n_split, -1);
  double exact_ic = 0;
  double ent_a = 0, ent_b = 0;
  for (int i = 0; i < a.n_split; ++i) ent_a += IcMatching(a.in[i], n - a.in[i], n);
  for (int j = 0; j < b.n_split; ++j) ent_b += IcMatching(b.in[j], n - b.in[j], n);
  p.entropy_bound = std::min(ent_a, ent_b) / n;
  for (int i = 0; i < a.n_split; ++i) {
    const int na = a.in[i], nA = n - na;
    for (int j = 0; j < b.n_split; ++j) {
      int ab = 0;
      for (int w = 0; w < a.nw; ++w) ab += __builtin_popcountll(a.row(i)[w] & b.row(j)[w]);
      const int nb = b.in[j];
      const int aB = na - ab, Ab = nb - ab, AB = nA - Ab;
      if ((!aB && !Ab) || (!ab && !AB)) {
        exact_ic += IcMatching(na, nA, n);
        a_match[i] = j;
        b_match[j] = i;
        break;
      }
    }
  }
  p.exact_part = exact_ic / n;
  p.rows.clear();
  p.cols.clear();
  for (int i = 0; i < a.n_split; ++i) if (a_match[i] < 0) p.rows.push_back(i);
  for (int j = 0; j < b.n_split; ++j) if (b_match[j] < 0) p.cols.push_back(j);
  p.d = std::max(p.rows.size(), p.cols.size());
  while (int(p.rows.size()) < p.d) p.rows.push_back(-1);
  while (int(p.cols.size()) < p.d) p.cols.push_back(-1);
  p.C.resize(p.d);
  for (int r = 0; r < p.d; ++r) {
    const int i = p.rows[r];
    for (int c = 0; c < p.d; ++c) {
      const int j = p.cols[c];
      if (i < 0 || j < 0) { p.C(r, c) = BIG; continue; }
      int ab = 0;
      for (int w = 0; w < a.nw; ++w) ab += __builtin_popcountll(a.row(i)[w] & b.row(j)[w]);
      p.C(r, c) = CellCost(a.in[i], b.in[j], ab, n);
    }
  }
}

static inline double ToScore(const Problem &p, i128 lb) {
  // Same arithmetic as mutual_clustering<true>(): (BIG * dim - cost) / BIG
  return double(cost(i128(BIG) * p.d - lb)) * (1.0 / double(BIG)) + p.exact_part;
}

// LB(v) = sum v + sum_i min_j (c_ij - v_j): valid lower bound on LAP cost.
static i128 LowerBound(const Problem &p, const std::vector<cost> &v,
                       std::vector<cost> *u_out = nullptr) {
  i128 lb = 0;
  for (int j = 0; j < p.d; ++j) lb += v[j];
  if (u_out) u_out->resize(p.d);
  for (int i = 0; i < p.d; ++i) {
    const cost *row = p.C.row(i);
    cost m = std::numeric_limits<cost>::max();
    for (int j = 0; j < p.d; ++j) m = std::min(m, row[j] - v[j]);
    lb += m;
    if (u_out) (*u_out)[i] = m;
  }
  return lb;
}

// From row duals u: v_j = min_i (c_ij - u_i).  LB = sum u + sum v.
static i128 LowerBoundU(const Problem &p, const std::vector<cost> &u,
                        std::vector<cost> *v_out = nullptr) {
  i128 lb = 0;
  for (int i = 0; i < p.d; ++i) lb += u[i];
  std::vector<cost> v(p.d, std::numeric_limits<cost>::max());
  for (int i = 0; i < p.d; ++i) {
    const cost *row = p.C.row(i);
    for (int j = 0; j < p.d; ++j) v[j] = std::min(v[j], row[j] - u[i]);
  }
  for (int j = 0; j < p.d; ++j) lb += v[j];
  if (v_out) *v_out = v;
  return lb;
}

// Alternating (block coordinate) dual ascent from v, k rounds; returns LB.
static i128 Alternate(const Problem &p, std::vector<cost> v, int k) {
  std::vector<cost> u;
  i128 lb = LowerBound(p, v, &u);
  for (int r = 0; r < k; ++r) {
    lb = std::max(lb, LowerBoundU(p, u, &v));
    lb = std::max(lb, LowerBound(p, v, &u));
  }
  return lb;
}

struct LapResult {
  cost opt = 0;
  std::vector<cost> v; // final column duals
  std::vector<int> rowsol;
  ana::LapWork work;
  uint64_t ticks = 0;
  // Augmentation-phase bound: augmentations done when LB first rejects
  int aug_total = 0, aug_needed = -1;
  double ub_rt = 0, ub_arr = 0, ub_col = 0;
  int ub_row_dummy = 0;
};

static LapResult Solve(Problem &p, const cost *v_init, bool track,
                       double min_score) {
  LapResult res;
  std::vector<lap_col> rowsol(p.d);
  std::vector<lap_row> colsol(p.d);
  LapScratch scratch;
  CostMatrix work_c = p.C; // lap may transpose
  if (!track) {
    const uint64_t t0 = __rdtsc();
    res.opt = ana::lap_hooked(p.d, work_c, rowsol, colsol, scratch, v_init,
                              [](int, const std::vector<cost> &, int, int) {},
                              res.work);
    res.ticks = __rdtsc() - t0;
  } else {
    std::vector<cost> vcopy;
    auto hook = [&](int phase, const std::vector<cost> &v, int done, int total) {
      vcopy.assign(v.begin(), v.begin() + p.d);
      if (phase == 0) res.ub_col = ToScore(p, LowerBound(p, vcopy));
      else if (phase == 1) res.ub_rt = ToScore(p, LowerBound(p, vcopy));
      else if (phase == 2) {
        res.ub_arr = ToScore(p, LowerBound(p, vcopy));
        res.aug_total = total;
        if (res.ub_arr < min_score) res.aug_needed = 0;
      } else if (phase == 3 && res.aug_needed < 0) {
        // Check at ~10% steps to limit analysis cost
        const int step = std::max(1, total / 10);
        if (done % step == 0 || done == total) {
          if (ToScore(p, LowerBound(p, vcopy)) < min_score) res.aug_needed = done;
        }
      }
    };
    res.opt = ana::lap_hooked(p.d, work_c, rowsol, colsol, scratch, v_init,
                              hook, res.work);
  }
  res.v.assign(scratch.v.begin(), scratch.v.begin() + p.d);
  res.rowsol.assign(rowsol.begin(), rowsol.end());
  return res;
}

// ---- Mapping splits between leaf sets ----
using Key = std::vector<uint64_t>;
struct KeyHash {
  size_t operator()(const Key &k) const {
    uint64_t h = 1469598103934665603ULL;
    for (auto w : k) { h ^= w; h *= 1099511628211ULL; h ^= h >> 29; }
    return h;
  }
};

// A set of splits on parent leaf indexing with associated dual values.
struct PriceSource {
  int nwp = 0; // words in parent indexing
  std::vector<Key> sets;
  std::vector<cost> val;
};

static Key ToParent(const SplitSet &s, int i, const std::vector<int> &map, int nwp) {
  Key k(nwp, 0);
  for (int t = 0; t < s.n_tip; ++t) {
    if ((s.row(i)[t / 64] >> (t % 64)) & 1) {
      const int pt = map[t];
      k[pt / 64] |= uint64_t(1) << (pt % 64);
    }
  }
  return k;
}

// Restrict to leaf set L, polarized so that reference leaf r is outside.
static Key Restrict(const Key &s, const Key &L, int r) {
  Key k(s.size());
  const bool flip = (s[r / 64] >> (r % 64)) & 1;
  for (size_t w = 0; w < s.size(); ++w) k[w] = (flip ? ~s[w] : s[w]) & L[w];
  return k;
}

enum Agg { AGG_MIN_V, AGG_MAX_V, AGG_MEAN_V };

// Map source prices onto the target problem's columns (if cols) or rows.
// Returns per-index value; `mapped` flags which received one.
static void MapPrices(const PriceSource &src, const SplitSet &target,
                      const std::vector<int> &target_idx,
                      const std::vector<int> &map, Agg agg,
                      std::vector<cost> &out, std::vector<char> &mapped,
                      int &n_multi) {
  const int nwp = src.nwp;
  Key L(nwp, 0);
  for (int t = 0; t < target.n_tip; ++t) L[map[t] / 64] |= uint64_t(1) << (map[t] % 64);
  const int r = map[0];
  std::unordered_map<Key, std::vector<int>, KeyHash> table;
  for (size_t q = 0; q < src.sets.size(); ++q) {
    table[Restrict(src.sets[q], L, r)].push_back(int(q));
  }
  const int d = target_idx.size();
  out.assign(d, 0);
  mapped.assign(d, 0);
  n_multi = 0;
  for (int c = 0; c < d; ++c) {
    const int j = target_idx[c];
    if (j < 0) continue;
    const Key k = Restrict(ToParent(target, j, map, nwp), L, r);
    auto it = table.find(k);
    if (it == table.end()) continue;
    const auto &qs = it->second;
    if (qs.size() > 1) ++n_multi;
    cost acc;
    if (agg == AGG_MIN_V) {
      acc = std::numeric_limits<cost>::max();
      for (int q : qs) acc = std::min(acc, src.val[q]);
    } else if (agg == AGG_MAX_V) {
      acc = std::numeric_limits<cost>::min();
      for (int q : qs) acc = std::max(acc, src.val[q]);
    } else {
      long double s = 0;
      for (int q : qs) s += src.val[q];
      acc = cost(s / qs.size());
    }
    out[c] = acc;
    mapped[c] = 1;
  }
}

// Fill unmapped column prices with the best feasible value given rows'
// minima over mapped columns: v_j = min_i (c_ij - u_i).
static bool CompleteV(const Problem &p, std::vector<cost> &v,
                      const std::vector<char> &mapped) {
  int nm = 0;
  for (char m : mapped) nm += m;
  if (nm == 0) return false;
  if (nm == p.d) return true;
  std::vector<cost> u(p.d);
  for (int i = 0; i < p.d; ++i) {
    cost m = std::numeric_limits<cost>::max();
    for (int j = 0; j < p.d; ++j) if (mapped[j]) m = std::min(m, p.C(i, j) - v[j]);
    u[i] = m;
  }
  for (int j = 0; j < p.d; ++j) {
    if (mapped[j]) continue;
    cost m = std::numeric_limits<cost>::max();
    for (int i = 0; i < p.d; ++i) m = std::min(m, p.C(i, j) - u[i]);
    v[j] = m;
  }
  return true;
}

static bool CompleteU(const Problem &p, std::vector<cost> &u,
                      const std::vector<char> &mapped) {
  int nm = 0;
  for (char m : mapped) nm += m;
  if (nm == 0) return false;
  if (nm == p.d) return true;
  std::vector<cost> v(p.d);
  for (int j = 0; j < p.d; ++j) {
    cost m = std::numeric_limits<cost>::max();
    for (int i = 0; i < p.d; ++i) if (mapped[i]) m = std::min(m, p.C(i, j) - u[i]);
    v[j] = m;
  }
  for (int i = 0; i < p.d; ++i) {
    if (mapped[i]) continue;
    cost m = std::numeric_limits<cost>::max();
    for (int j = 0; j < p.d; ++j) m = std::min(m, p.C(i, j) - v[j]);
    u[i] = m;
  }
  return true;
}

// Size-only bound: max-weight matching where each weight is the largest
// MCI attainable by two splits of these sizes (nested configurations).
static double SizeBound(const Problem &p, const SplitSet &a, const SplitSet &b) {
  Problem q;
  q.n_tip = p.n_tip; q.d = p.d; q.exact_part = p.exact_part;
  q.C.resize(p.d);
  const int n = p.n_tip;
  for (int r = 0; r < p.d; ++r) {
    for (int c = 0; c < p.d; ++c) {
      const int i = p.rows[r], j = p.cols[c];
      if (i < 0 || j < 0) { q.C(r, c) = BIG; continue; }
      const int na = a.in[i], nb = b.in[j];
      const cost c1 = CellCost(na, nb, std::min(na, nb), n);
      const cost c2 = CellCost(na, nb, std::max(0, na + nb - n), n);
      q.C(r, c) = std::min(c1, c2);
    }
  }
  LapResult res = Solve(q, nullptr, false, 0);
  return ToScore(q, res.opt);
}

// Progressive priced bound: rows filled in `order`; unfilled rows use the
// provable bound c_ij >= BIG - floor(BIG * min(h_i, h_j)) - 1.  Returns the
// fraction of rows that must be filled before the bound rejects (2 = never).
static double ProgressiveFraction(const Problem &p, const SplitSet &a,
                                  const SplitSet &b, const std::vector<cost> &v,
                                  double min_score, int order_kind) {
  const int n = p.n_tip, d = p.d;
  std::vector<double> hr(d, 0), hc(d, 0);
  for (int r = 0; r < d; ++r) if (p.rows[r] >= 0) hr[r] = IcMatching(a.in[p.rows[r]], n - a.in[p.rows[r]], n) / n;
  for (int c = 0; c < d; ++c) if (p.cols[c] >= 0) hc[c] = IcMatching(b.in[p.cols[c]], n - b.in[p.cols[c]], n) / n;
  std::vector<i128> exact_u(d), loose_u(d);
  i128 lb = 0;
  for (int j = 0; j < d; ++j) lb += v[j];
  for (int i = 0; i < d; ++i) {
    cost mu = std::numeric_limits<cost>::max(), ml = std::numeric_limits<cost>::max();
    for (int j = 0; j < d; ++j) {
      mu = std::min(mu, p.C(i, j) - v[j]);
      const cost lower = BIG - cost(std::floor(double(BIG) * std::min(hr[i], hc[j]))) - 1;
      ml = std::min(ml, std::max(cost(0), lower) - v[j]);
    }
    exact_u[i] = mu;
    loose_u[i] = ml;
    lb += ml;
  }
  std::vector<int> ord(d);
  for (int i = 0; i < d; ++i) ord[i] = i;
  if (order_kind == 1) {
    // Oracle-free heuristic: largest provable slack first
    std::sort(ord.begin(), ord.end(), [&](int x, int y) {
      return hr[x] > hr[y];
    });
  } else if (order_kind == 2) {
    // Oracle: largest actual tightening first (best possible order)
    std::sort(ord.begin(), ord.end(), [&](int x, int y) {
      return exact_u[x] - loose_u[x] > exact_u[y] - loose_u[y];
    });
  }
  if (ToScore(p, lb) < min_score) return 0;
  for (int t = 0; t < d; ++t) {
    const int i = ord[t];
    lb += exact_u[i] - loose_u[i];
    if (ToScore(p, lb) < min_score) return double(t + 1) / d;
  }
  return 2;
}

struct Child {
  int cand, kind, n_tip;
  std::vector<int> map;
  double min_score, exact, bounded, info;
  SplitSet a, b;
};

int main(int argc, char **argv) {
  if (argc < 3) { fprintf(stderr, "usage\n"); return 1; }
  InitLg2();
  FILE *f = fopen(argv[1], "rb");
  FILE *out = fopen(argv[2], "w");
  const bool do_size = argc > 3 && std::string(argv[3]) == "size";
  fprintf(out, "parent,cand,kind,n_tip,d,exact,lap_exact,min_score,"
               "ub_entropy,ub_row,ub_cr,ub_rc,ub_rt,ub_arr,"
               "ub_alt1,ub_alt3,ub_alt10,"
               "ub_pv_minv,ub_pv_maxv,ub_pv_mean,ub_pu,ub_pv_alt1,ub_pv_alt3,"
               "ub_pv_warmarr,ub_prev,ub_size,frac_mapped,n_multi,"
               "aug_total,aug_needed,"
               "cold_free_init,cold_free_arr,cold_arr,cold_aug,cold_ticks,"
               "warm_free_init,warm_free_arr,warm_arr,warm_aug,warm_ticks,"
               "warm_ok,prev_ticks,parent_n_tip,prog_nat,prog_h,prog_oracle\n");
  int parent_id = -1;
  int parent_n_tip = 0;
  PriceSource parent_v, parent_u;
  PriceSource prev_v;
  bool have_prev = false;
  int prev_parent = -2;
  long n_records = 0, n_mismatch = 0;
  int tag;
  while (fread(&tag, sizeof(int), 1, f) == 1) {
    if (tag == 1) {
      int n_tip;
      double info;
      fread(&n_tip, sizeof(int), 1, f);
      fread(&info, sizeof(double), 1, f);
      SplitSet A, B;
      ReadSplits(f, n_tip, A);
      ReadSplits(f, n_tip, B);
      ++parent_id;
      parent_n_tip = n_tip;
      have_prev = false;
      Problem p;
      Build(A, B, p);
      const int nwp = (n_tip + 63) / 64;
      std::vector<int> ident(n_tip);
      for (int t = 0; t < n_tip; ++t) ident[t] = t;
      parent_v = PriceSource{nwp, {}, {}};
      parent_u = PriceSource{nwp, {}, {}};
      if (p.d > 1) {
        LapResult res = Solve(p, nullptr, false, 0);
        // Check optimality of the final duals
        std::vector<cost> u;
        const i128 lb = LowerBound(p, res.v, &u);
        if (lb != i128(res.opt)) {
          fprintf(stderr, "parent %d: dual gap %g\n", parent_id, double(i128(res.opt) - lb));
        }
        {
          cost mn = *std::min_element(res.v.begin(), res.v.end());
          cost mx = *std::max_element(res.v.begin(), res.v.end());
          static double worst = 0;
          const double spread = double(mx - mn) / double(BIG);
          if (spread > worst) { worst = spread; fprintf(stderr, "parent %d d=%d: dual spread %.3f BIG\n", parent_id, p.d, spread); }
        }
        for (int c = 0; c < p.d; ++c) {
          if (p.cols[c] < 0) continue;
          parent_v.sets.push_back(ToParent(B, p.cols[c], ident, nwp));
          parent_v.val.push_back(res.v[c]);
        }
        for (int r = 0; r < p.d; ++r) {
          if (p.rows[r] < 0) continue;
          parent_u.sets.push_back(ToParent(A, p.rows[r], ident, nwp));
          parent_u.val.push_back(u[r]);
        }
      }
      continue;
    }
    if (tag != 2) { fprintf(stderr, "bad tag %d\n", tag); return 1; }
    Child ch;
    int hdr[3];
    fread(hdr, sizeof(int), 3, f);
    ch.cand = hdr[0]; ch.kind = hdr[1]; ch.n_tip = hdr[2];
    ch.map.resize(ch.n_tip);
    fread(ch.map.data(), sizeof(int), ch.n_tip, f);
    double vals[4];
    fread(vals, sizeof(double), 4, f);
    ch.min_score = vals[0]; ch.exact = vals[1]; ch.bounded = vals[2]; ch.info = vals[3];
    ReadSplits(f, ch.n_tip, ch.a);
    ReadSplits(f, ch.n_tip, ch.b);
    ++n_records;

    Problem p;
    Build(ch.a, ch.b, p);
    if (p.d == 0 || ch.a.n_split == 0 || ch.b.n_split == 0 ||
        int(p.rows.size()) == 0) {
      continue;
    }
    // All-exact case
    bool all_exact = true;
    for (int r = 0; r < p.d; ++r) if (p.rows[r] >= 0 && p.cols[r] >= 0) all_exact = false;
    if (p.d == 0 || all_exact) continue;

    LapResult cold = Solve(p, nullptr, true, ch.min_score);
    const double lap_exact = ToScore(p, cold.opt);
    if (std::fabs(lap_exact - ch.exact) > 1e-9) ++n_mismatch;

    // Row and simple bounds
    std::vector<cost> zero(p.d, 0), u0;
    const double ub_row = ToScore(p, LowerBound(p, zero, &u0));
    std::vector<cost> v_rc;
    const double ub_rc = ToScore(p, LowerBoundU(p, u0, &v_rc));
    // column-first then alternating
    std::vector<cost> vcol(p.d, std::numeric_limits<cost>::max());
    for (int i = 0; i < p.d; ++i) for (int j = 0; j < p.d; ++j) vcol[j] = std::min(vcol[j], p.C(i, j));
    const double ub_alt1 = ToScore(p, std::max(Alternate(p, vcol, 1), LowerBoundU(p, u0)));
    const double ub_alt3 = ToScore(p, Alternate(p, vcol, 3));
    const double ub_alt10 = ToScore(p, Alternate(p, vcol, 10));

    // Parent prices
    double ub_pv[3] = {NAN, NAN, NAN}, ub_pu = NAN, ub_pv_alt1 = NAN, ub_pv_alt3 = NAN;
    double ub_warmarr = NAN, frac_mapped = NAN;
    int n_multi = 0;
    LapResult warm;
    bool warm_ok = false;
    std::vector<cost> v_best;
    if (!parent_v.sets.empty()) {
      for (int a = 0; a < 3; ++a) {
        std::vector<cost> v;
        std::vector<char> mapped;
        int nm;
        MapPrices(parent_v, ch.b, p.cols, ch.map, Agg(a), v, mapped, nm);
        if (a == 0) {
          int cnt = 0;
          for (char m : mapped) cnt += m;
          frac_mapped = double(cnt) / p.d;
          n_multi = nm;
        }
        if (CompleteV(p, v, mapped)) {
          ub_pv[a] = ToScore(p, LowerBound(p, v));
          if (a == 0) v_best = v;
        }
      }
      {
        std::vector<cost> u;
        std::vector<char> mapped;
        int nm;
        MapPrices(parent_u, ch.a, p.rows, ch.map, AGG_MAX_V, u, mapped, nm);
        if (CompleteU(p, u, mapped)) ub_pu = ToScore(p, LowerBoundU(p, u));
      }
      if (!v_best.empty()) {
        ub_pv_alt1 = ToScore(p, Alternate(p, v_best, 1));
        ub_pv_alt3 = ToScore(p, Alternate(p, v_best, 3));
        // Warm-started JV from parent prices; bound after ARR
        double ub_after_arr = NAN;
        ana::LapWork w;
        LapResult wr;
        {
          std::vector<lap_col> rowsol(p.d);
          std::vector<lap_row> colsol(p.d);
          LapScratch scratch;
          CostMatrix work_c = p.C;
          std::vector<cost> vv;
          auto hook = [&](int phase, const std::vector<cost> &v, int, int) {
            if (phase == 2) {
              vv.assign(v.begin(), v.begin() + p.d);
              ub_after_arr = ToScore(p, LowerBound(p, vv));
            }
          };
          wr.opt = ana::lap_hooked(p.d, work_c, rowsol, colsol, scratch,
                                   v_best.data(), hook, wr.work);
        }
        ub_warmarr = ub_after_arr;
        // Timed warm run without hooks
        warm = Solve(p, v_best.data(), false, 0);
        warm_ok = warm.opt == cold.opt;
        if (!warm_ok) {
          const i128 lbw = LowerBound(p, warm.v), lbc = LowerBound(p, cold.v);
          fprintf(stderr, "mismatch d=%d warm-cold %lld; warm gap %lld; cold gap %lld (units; 1 bit = %lld)\n",
                  p.d, (long long)(warm.opt - cold.opt), (long long)(warm.opt - lbw),
                  (long long)(cold.opt - lbc), (long long)BIG);
        }
      }
    }
    // Previous call's optimal prices (same parent)
    double ub_prev = NAN;
    unsigned long long prev_ticks = 0;
    if (have_prev && prev_parent == parent_id) {
      std::vector<cost> v;
      std::vector<char> mapped;
      int nm;
      MapPrices(prev_v, ch.b, p.cols, ch.map, AGG_MIN_V, v, mapped, nm);
      if (CompleteV(p, v, mapped)) {
        ub_prev = ToScore(p, LowerBound(p, v));
        LapResult pr = Solve(p, v.data(), false, 0);
        prev_ticks = pr.ticks;
      }
    }
    // Remember this call's optimal duals
    {
      const int nwp = parent_v.nwp ? parent_v.nwp : (ch.n_tip + 63) / 64;
      prev_v = PriceSource{nwp, {}, {}};
      for (int c = 0; c < p.d; ++c) {
        if (p.cols[c] < 0) continue;
        prev_v.sets.push_back(ToParent(ch.b, p.cols[c], ch.map, nwp));
        prev_v.val.push_back(cold.v[c]);
      }
      have_prev = true;
      prev_parent = parent_id;
    }
    const double ub_size = do_size ? SizeBound(p, ch.a, ch.b) : NAN;
    double prog[3] = {NAN, NAN, NAN};
    if (!v_best.empty()) {
      for (int k = 0; k < 3; ++k) prog[k] = ProgressiveFraction(p, ch.a, ch.b, v_best, ch.min_score, k);
    }
    // Timed cold run without hooks
    LapResult cold_t = Solve(p, nullptr, false, 0);
    {
      const int kk = (cold.ub_row_dummy, 0);
      (void)kk;
      g_init[0] += cold_t.work.t_init; g_arr[0] += cold_t.work.t_arr; g_aug[0] += cold_t.work.t_aug; ++g_n[0];
      g_init[1] += warm.work.t_init; g_arr[1] += warm.work.t_arr; g_aug[1] += warm.work.t_aug; ++g_n[1];
    }

    fprintf(out, "%d,%d,%d,%d,%d,%.12g,%.12g,%.12g,"
                 "%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,"
                 "%.12g,%.12g,%.12g,"
                 "%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,"
                 "%.12g,%.12g,%.12g,%.4g,%d,"
                 "%d,%d,"
                 "%lld,%lld,%lld,%lld,%llu,"
                 "%lld,%lld,%lld,%lld,%llu,%d,%llu,%d,%.4g,%.4g,%.4g\n",
            parent_id, ch.cand, ch.kind, ch.n_tip, p.d, ch.exact, lap_exact,
            ch.min_score, p.entropy_bound, ub_row, cold.ub_col, ub_rc,
            cold.ub_rt, cold.ub_arr, ub_alt1, ub_alt3, ub_alt10,
            ub_pv[0], ub_pv[1], ub_pv[2], ub_pu, ub_pv_alt1, ub_pv_alt3,
            ub_warmarr, ub_prev, ub_size, frac_mapped, n_multi,
            cold.aug_total, cold.aug_needed,
            cold_t.work.n_free_after_init, cold_t.work.n_free_after_arr,
            cold_t.work.arr_steps, cold_t.work.aug_inner + cold_t.work.aug_scan,
            (unsigned long long)cold_t.ticks,
            warm.work.n_free_after_init, warm.work.n_free_after_arr,
            warm.work.arr_steps, warm.work.aug_inner + warm.work.aug_scan,
            (unsigned long long)warm.ticks, int(warm_ok), prev_ticks, parent_n_tip,
            prog[0], prog[1], prog[2]);
  }
  fprintf(stderr, "records %ld, exact mismatches %ld\n", n_records, n_mismatch);
  for (int k = 0; k < 2; ++k) {
    const double tot = double(g_init[k] + g_arr[k] + g_aug[k]);
    fprintf(stderr, "%s LAP: mean %.1f us; init %.1f%%, ARR %.1f%%, augment %.1f%%\n", k ? "warm" : "cold",
            tot / g_n[k] / 2100.0, 100 * g_init[k] / tot, 100 * g_arr[k] / tot, 100 * g_aug[k] / tot);
  }
  fclose(out);
  fclose(f);
}
