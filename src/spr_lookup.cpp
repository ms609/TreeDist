#include <Rcpp.h>
#include <algorithm>
#include <array>
#include "spr/lookup_table_7.h"

using namespace Rcpp;

inline int count_bits(uint8_t n) { return __builtin_popcount(n); }

// Standard normalization: ensure split has < 4 tips and bit 0 is 0
inline uint8_t norm(uint8_t s) {
  if (count_bits(s) > 3) s = ~s & 0x7F;
  return s;
}

template <size_t N>
double search_table(const std::array<SPRScore, N>& table, uint32_t key) {
  auto it = std::lower_bound(table.begin(), table.end(), key, 
                             [](const SPRScore& a, uint32_t b) {
                               return a.key < b;
                             });
  if (it != table.end() && it->key == key) {
    return static_cast<double>(it->score);
  }
  return NA_REAL;
}

void get_pec_map(const std::array<uint8_t, 4>& s1, int* map) {
  uint8_t t[2], p[2];
  int ti = 0, pi = 0;
  
  for (uint8_t s : s1) {
    uint8_t n = norm(s);
    if (count_bits(n) == 3) t[ti++] = n;
    else p[pi++] = n;
  }
  
  uint8_t midpoint = t[0] ^ t[1];
  // R: if (!TipsInSplits(xor(p1, t1)) %in% c(1, 6)) { pairs <- pairs[2:1] }
  if (count_bits(p[0] ^ t[0]) != 1) std::swap(p[0], p[1]);
  
  int midTip = __builtin_ctz(midpoint);
  uint8_t trio1Tip = norm(t[0] ^ p[0]);
  uint8_t trio2Tip = norm(t[1] ^ p[1]);
  uint8_t duo1Tips = p[0];
  uint8_t duo2Tips = p[1];
  
  // Filling the map based on canonOrder:
  // c(midTip, which(trio1Tip), which(duo1Tips), which(trio2Tip), which(duo2Tips))
  map[midTip] = 0;
  map[__builtin_ctz(trio1Tip)] = 1;
  
  // Duo1 has two tips
  int d1a = __builtin_ctz(duo1Tips);
  int d1b = __builtin_ctz(duo1Tips ^ (1 << d1a));
  map[d1a] = 2; map[d1b] = 3;
  
  map[__builtin_ctz(trio2Tip)] = 4;
  
  // Duo2 has two tips
  int d2a = __builtin_ctz(duo2Tips);
  int d2b = __builtin_ctz(duo2Tips ^ (1 << d2a));
  map[d2a] = 5; map[d2b] = 6;
}

// [[Rcpp::export]]
double lookup_from_table(RawVector sp1, RawVector sp2) {
  std::array<uint8_t, 4> s1, s2;
  int trios = 0;
  for (int i = 0; i < 4; ++i) {
    s1[i] = sp1[i];
    s2[i] = sp2[i];
    if (count_bits(norm(s1[i])) == 3) trios++;
  }
  
  int map[7];
  bool is_pec = (trios == 2);
  
  if (is_pec) {
    get_pec_map(s1, map);
  } else {
    // ... Similar logic for get_bal_map ...
  }
  
  std::array<int, 4> remapped_s2;
  for (int i = 0; i < 4; ++i) {
    uint8_t res = 0;
    for (int bit = 0; bit < 7; ++bit) {
      if (s2[i] & (1 << bit)) res |= (1 << map[bit]);
    }
    remapped_s2[i] = (int)norm(res);
  }
  std::sort(remapped_s2.begin(), remapped_s2.end());
  
  uint32_t key = (static_cast<uint32_t>(remapped_s2[0] - 6) << 20) |
    (static_cast<uint32_t>(remapped_s2[1] - 14) << 13) |
    (static_cast<uint32_t>(remapped_s2[2] - 30) << 6)  |
    (static_cast<uint32_t>(remapped_s2[3] - 62));
  
  if (is_pec) {
    return search_table(PEC_LOOKUP, key);
  } else {
    return search_table(BAL_LOOKUP, key);
  }
}