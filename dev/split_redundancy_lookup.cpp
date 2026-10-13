// Fill an n x n cost block either by computing MCI cell cost from split bits,
// or by gathering from a precomputed U x U table indexed by unique-split ids.
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <random>
#include <vector>
int main(int argc, char** argv) {
  const int nTip = atoi(argv[1]), U = atoi(argv[2]), n = nTip - 3, reps = 2000;
  const int bins = (nTip + 63) / 64;
  std::mt19937_64 rng(1);
  std::vector<uint64_t> bits((size_t)U * bins);
  std::vector<int> inSplit(U);
  for (int u = 0; u < U; ++u) { int c = 0;
    for (int b = 0; b < bins; ++b) { bits[(size_t)u*bins+b] = rng(); c += __builtin_popcountll(bits[(size_t)u*bins+b]); }
    inSplit[u] = c; }
  std::vector<double> lg2(nTip + 1); for (int i = 1; i <= nTip; ++i) lg2[i] = std::log2(i);
  std::vector<int64_t> table((size_t)U * U);
  for (auto& x : table) x = rng() >> 8;
  std::vector<int64_t> out((size_t)n * n);
  std::uniform_int_distribution<int> pick(0, U - 1);
  std::vector<int> a(n), b(n);
  double tC = 0, tL = 0; volatile int64_t sink = 0;
  for (int r = 0; r < reps; ++r) {
    for (int i = 0; i < n; ++i) { a[i] = pick(rng); b[i] = pick(rng); }
    auto t0 = std::chrono::steady_clock::now();
    for (int i = 0; i < n; ++i) { const uint64_t* ar = &bits[(size_t)a[i]*bins]; int na = inSplit[a[i]], nA = nTip - na;
      for (int j = 0; j < n; ++j) { const uint64_t* br = &bits[(size_t)b[j]*bins]; int ab = 0;
        for (int k = 0; k < bins; ++k) ab += __builtin_popcountll(ar[k] & br[k]);
        int nb = inSplit[b[j]], nB = nTip - nb, aB = na - ab, Ab = nb - ab, AB = nA - Ab;
        double s = ab * (lg2[ab] - lg2[na] - lg2[nb]) + aB * (lg2[aB] - lg2[na] - lg2[nB]) +
                   Ab * (lg2[Ab] - lg2[nA] - lg2[nb]) + AB * (lg2[AB] - lg2[nA] - lg2[nB]);
        out[(size_t)i*n+j] = (int64_t)(s * 1e6); } }
    auto t1 = std::chrono::steady_clock::now();
    sink += out[r % (n*n)];
    for (int i = 0; i < n; ++i) { const int64_t* row = &table[(size_t)a[i]*U];
      for (int j = 0; j < n; ++j) out[(size_t)i*n+j] = row[b[j]]; }
    auto t2 = std::chrono::steady_clock::now();
    sink += out[r % (n*n)];
    tC += std::chrono::duration<double>(t1 - t0).count(); tL += std::chrono::duration<double>(t2 - t1).count();
  }
  printf("nTip=%d U=%d tableMB=%.0f  compute=%.2f ns/cell  lookup=%.2f ns/cell\n", nTip, U, U*(double)U*8/1e6,
         tC/reps/n/n*1e9, tL/reps/n/n*1e9);
}
