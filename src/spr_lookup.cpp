#include <Rcpp.h>
#include <algorithm>
#include <array>
#include "spr/lookup_table_7.h"

using namespace Rcpp;

inline int count_bits(uint8_t n) { return __builtin_popcount(n); }

// Polarize on Tip 7 (bit 6): Ensure bit 6 is always 0
// This matches R: !PolarizeSplits(splits, 7)
inline uint8_t norm7(uint8_t s) {
  if (s & 0x40) { // If tip 7 (bit 6) is set
    return ~s & 0x7F; // Flip bits and mask to 7 bits
  }
  return s & 0x7F;
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

// Logic to identify the map from sp1 to canonOrder
void get_pec_map(const std::array<uint8_t, 4>& s1, int* map) {
  uint8_t t[2], p[2];
  int ti = 0, pi = 0;
  
  for (uint8_t s : s1) {
    uint8_t n = norm7(s);
    if (count_bits(n) == 3) t[ti++] = n;
    else p[pi++] = n;
  }
  
  uint8_t midpoint = t[0] ^ t[1];
  // Align pair1 with trio1
  if (count_bits(p[0] ^ t[0]) != 1) std::swap(p[0], p[1]);
  
  int midTip = __builtin_ctz(midpoint);
  uint8_t trio1Tip = norm7(t[0] ^ p[0]);
  uint8_t trio2Tip = norm7(t[1] ^ p[1]);
  uint8_t duo1Tips = p[0];
  uint8_t duo2Tips = p[1];
  
  map[midTip] = 0;
  map[__builtin_ctz(trio1Tip)] = 1;
  
  int d1a = __builtin_ctz(duo1Tips);
  map[d1a] = 2; map[__builtin_ctz(duo1Tips ^ (1 << d1a))] = 3;
  
  map[__builtin_ctz(trio2Tip)] = 4;
  
  int d2a = __builtin_ctz(duo2Tips);
  map[d2a] = 5; map[__builtin_ctz(duo2Tips ^ (1 << d2a))] = 6;
}

void get_bal_map(const std::array<uint8_t, 4>& s1, int* map) {
  int tiss[4];
  for (int i = 0; i < 4; ++i) tiss[i] = count_bits(norm7(s1[i]));
  
  int firstTrioIdx = 0;
  for (int i = 1; i < 4; ++i) if (tiss[i] > tiss[firstTrioIdx]) firstTrioIdx = i;
  
  uint8_t firstTrioSp = norm7(s1[firstTrioIdx]);
  
  uint8_t soloSp = 0;
  int trioPairIdx = -1;
  for (int i = 0; i < 4; ++i) {
    if (i == firstTrioIdx) continue;
    uint8_t curSp = norm7(s1[i]);
    uint8_t testSolo = firstTrioSp ^ curSp;
    if (count_bits(testSolo) == 6 || count_bits(testSolo) == 1) {
      soloSp = norm7(testSolo);
      trioPairIdx = i;
      break;
    }
  }
  
  int other[2], o_count = 0;
  for (int i = 0; i < 4; ++i) {
    if (i != firstTrioIdx && i != trioPairIdx) other[o_count++] = i;
  }
  
  uint8_t singleton = norm7(~soloSp & 0x7F);
  int singleTip = __builtin_ctz(singleton);
  uint8_t trioPairTip = norm7(s1[trioPairIdx]);
  uint8_t osp1 = norm7(s1[other[0]]);
  uint8_t osp2 = norm7(s1[other[1]]);
  
  map[singleTip] = 0;
  int tp1 = __builtin_ctz(trioPairTip);
  map[tp1] = 1; map[__builtin_ctz(trioPairTip ^ (1 << tp1))] = 2;
  int o1a = __builtin_ctz(osp1);
  map[o1a] = 3; map[__builtin_ctz(osp1 ^ (1 << o1a))] = 4;
  int o2a = __builtin_ctz(osp2);
  map[o2a] = 5; map[__builtin_ctz(osp2 ^ (1 << o2a))] = 6;
}

// [[Rcpp::export]]
double lookup_from_table(RawVector sp1, RawVector sp2) {
  if (sp1.size() < 4 || sp2.size() < 4) return NA_REAL;
  
  std::array<uint8_t, 4> s1, s2;
  int trios = 0;
  for (int i = 0; i < 4; ++i) {
    s1[i] = sp1[i];
    s2[i] = sp2[i];
    if (count_bits(norm7(s1[i])) == 3) trios++;
  }
  
  int map[7];
  bool is_pec = (trios == 2);
  
  if (is_pec) get_pec_map(s1, map);
  else get_bal_map(s1, map);
  
  std::array<int, 4> remapped_s2;
  for (int i = 0; i < 4; ++i) {
    uint8_t res = 0;
    for (int bit = 0; bit < 7; ++bit) {
      if (s2[i] & (1 << bit)) res |= (1 << map[bit]);
    }
    remapped_s2[i] = (int)norm7(res);
  }
  std::sort(remapped_s2.begin(), remapped_s2.end());
  
  // Check if within bounds for the ultra-compact pack
  if (remapped_s2[0] < 3 || remapped_s2[1] < 7 || 
      remapped_s2[2] < 15 || remapped_s2[3] < 33) return NA_REAL;
  
  uint32_t key = (static_cast<uint32_t>(remapped_s2[0] - 3) << 18) |
    (static_cast<uint32_t>(remapped_s2[1] - 7) << 12) |
    (static_cast<uint32_t>(remapped_s2[2] - 15) << 6) |
    (static_cast<uint32_t>(remapped_s2[3] - 33));
  
  return is_pec ? search_table(PEC_LOOKUP, key) : search_table(BAL_LOOKUP, key);
}
