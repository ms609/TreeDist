#include <Rcpp.h>
#include <algorithm>
#include <array>
#include <cstdint>
#include "spr/lookup.h"
#include "spr/canonize7.h"
#include "spr/canonize8.h"
#include "spr/canonize9.h"
#include "spr/lookup_table_7.h"
#include "spr/lookup_table_8.h"
#include "TreeTools/assert.h"


inline Shape7 detect_shape7(const SplitSet7& sp) {
  int trio_count = 0;
  for (auto s : sp) {
    const int t = popcount7(s);
    if (t == 3) ++trio_count;
  }
  return (trio_count == 2) ? Shape7::Pectinate : Shape7::Balanced;
}

inline Shape8 detect_shape8(const SplitSet8& sp) {
  int n4 = 0, n3 = 0, n2 = 0;
  for (auto s : sp) {
    const int t = tips_in_smallest8(s);
    if (t == 4) ++n4;
    else if (t == 3) ++n3;
    else if (t == 2) ++n2;
  }
  
  if (n4 == 1 && n3 == 2) return Shape8::Pectinate;
  if (n4 == 1 && n3 == 1) return Shape8::Mix;
  if (n3 == 2)            return Shape8::Mid;
  return Shape8::Balanced;
}

inline Shape9 detect_shape9(const SplitSet9& sp) {
  int n4 = 0, n3 = 0, n2 = 0;
  
  for (auto s : sp) {
    const int t = tips_in_smallest9(s);
    if      (t == 4) ++n4;
    else if (t == 3) ++n3;
    else if (t == 2) ++n2;
  }
  
  if (n4 == 2 && n3 == 2) return Shape9::s0;
  if (n4 == 2 && n3 == 1) return Shape9::s1;
  if (n4 == 1 && n3 == 2) return Shape9::s2;
  if (n4 == 1 && n3 == 1) return Shape9::s3;
  if (n3 == 3)            return Shape9::s4;
  return Shape9::s5;
}

// Reorder the leaves of s to match the sequence given by p
Split7 permute_split7(Split7 s, const Perm7& p) {
  Split7 out = 0;
  for (int i = 0; i < 7; ++i) {
    if (s & (1 << p[i])) {
      out |= (1 << i);
    }
  }
  return out;
}

Split8 permute_split8(Split8 s, const Perm8& p) {
  Split8 out = 0;
  for (int i = 0; i < 8; ++i) {
    if (s & (1 << p[i])) {
      out |= (1 << i);
    }
  }
  return out;
}

Split9 permute_split9(Split9 s, const Perm9& p) {
  Split9 out = 0;
  for (int i = 0; i < 9; ++i) {
    if (s & (1 << p[i])) {
      out |= (1 << i);
    }
  }
  return out;
}

inline Split7 polarize7(Split7 s) {
  if (s & (1 << 6)) s ^= MASK7;
  return s;
}

inline Split8 polarize8(Split8 s) {
  if (s & (1 << 7)) s ^= MASK8;
  return s;
}

inline Split9 polarize9(Split9 s) {
  if (s & (1 << 8)) s ^= MASK9;
  return s;
}

inline uint32_t BitPack7(const std::array<int,4>& v) {
  return ((v[0] - 3)  << 18) |
    ((v[1] - 7)  << 12) |
    ((v[2] - 15) << 6) |
    ( v[3] - 33);
}

inline uint64_t BitPack8(const std::array<int,5>& v) {
  return
  ((uint64_t)(v[0] - 3) << 27) |
    ((uint64_t)(v[1] - 7) << 20) |
    ((uint64_t)(v[2] - 15) << 13) |
    ((uint64_t)(v[3] - 31) << 6)  |
    (uint64_t)(v[4] - 65);
}

inline uint64_t BitPack9(const std::array<int,6>& v) {
  Rcpp::stop("We've not populated this yet!");
  return
  ((uint64_t)v[0] << 27) |
    ((uint64_t)v[1] << 20) |
    ((uint64_t)v[2] << 13) |
    ((uint64_t)v[3] << 6)  |
    ((uint64_t)v[4] << 6)  |   // intentional duplication
    (uint64_t)v[5];
}

template <size_t N>
int lookup7(uint32_t key, const std::array<SPRScore, N>& table) {
  auto it = std::lower_bound(
    table.begin(), table.end(), key,
    [](const SPRScore& a, uint32_t k) { return a.key < k; }
  );
  return (it != table.end() && it->key == key) ? it->score : -1;
}

template <size_t N>
int lookup8(uint64_t key, const std::array<SPRScore64, N>& table) {
  auto it = std::lower_bound(
    table.begin(), table.end(), key,
    [](const SPRScore64& a, uint64_t k) { return a.key < k; }
  );
  return (it != table.end() && it->key == key) ? it->score : -1;
}

int lookup6(const SplitSet6& sp1_raw, const SplitSet6& sp2_raw) {
  
  SplitSet6 sp1 = sp1_raw;
  SplitSet6 sp2 = sp2_raw;
  
  for (auto& s : sp1) s = smaller_split6(s);
  for (auto& s : sp2) s = smaller_split6(s);
  
  std::array<bool,3> pairs1{}, pairs2{};
  
  for (int i = 0; i < 3; ++i) {
    pairs1[i] = (popcount6(sp1[i]) == 2);
    pairs2[i] = (popcount6(sp2[i]) == 2);
  }
  
  /* All cherries shortcut */
  if ((pairs1[0] && pairs1[1] && pairs1[2]) ||
  (pairs2[0] && pairs2[1] && pairs2[2])) {
    return 2;
  }
  
  Split6 duo1[2], duo2[2];
  Split6 trio1 = 0, trio2 = 0;
  
  int d = 0;
  for (int i = 0; i < 3; ++i) {
    if (pairs1[i]) duo1[d++] = sp1[i];
    else trio1 = sp1[i];
  }
  
  d = 0;
  for (int i = 0; i < 3; ++i) {
    if (pairs2[i]) duo2[d++] = sp2[i];
    else trio2 = sp2[i];
  }
  
  Split6 middle1a = overlapper6(duo1[0], trio1);
  Split6 middle1b = overlapper6(duo1[1], trio1);
  
  Split6 middle2a = overlapper6(duo2[0], trio2);
  Split6 middle2b = overlapper6(duo2[1], trio2);
  
  Split6 inMiddle = (middle1a | middle1b) & (middle2a | middle2b);
  
  if (popcount6(inMiddle) == 1) {
    Split6 La = (middle1a & inMiddle) ? middle1a : middle1b;
    Split6 Lbc1 = (middle1a & inMiddle) ? duo1[0] : duo1[1];
    
    if (Lbc1 & La) Lbc1 ^= MASK6;
    if (Lbc1 & (middle2a | middle2b)) {
      
      Split6 Lbc2 = (middle2a & La) ? duo2[0] : duo2[1];
      
      if (Lbc2 & La) Lbc2 ^= MASK6;
      
      if ((Lbc2 & Lbc1) == 0) return 1;
    }
  }
  
  return 2;
}


int lookup7(const SplitSet7& sp1, const SplitSet7& sp2) {
  Shape7 shape = detect_shape7(sp1);
  
  CanonicalInfo7 canon =
    (shape == Shape7::Pectinate)
    ? canonical_pectinate(sp1)
      : canonical_balanced(sp1);
  
  std::array<int,4> packed{};
  for (int i = 0; i < 4; ++i) {
    Split7 s = permute_split7(sp2[i], canon.perm);
    s = polarize7(s);
    packed[i] = s;
  }
  
  std::sort(packed.begin(), packed.end());
  
  uint32_t key = BitPack7(packed);
  return (shape == Shape7::Pectinate)
    ? lookup7(key, PEC_LOOKUP7)
      : lookup7(key, BAL_LOOKUP7);
}

int lookup8(const SplitSet8& sp1, const SplitSet8& sp2) {
  Shape8 shape = detect_shape8(sp1);
  
  CanonicalInfo8 canon =
    (shape == Shape8::Pectinate) ? canonical_pectinate8(sp1) :
    (shape == Shape8::Mix)       ? canonical_mix8(sp1) :
    (shape == Shape8::Mid)       ? canonical_mid8(sp1) :
    canonical_balanced8(sp1);
  
  std::array<int,5> packed{};
  for (int i = 0; i < 5; ++i) {
    Split8 s = permute_split8(sp2[i], canon.perm);
    s = polarize8(s);
    packed[i] = s;
  }
  
  std::sort(packed.begin(), packed.end());
  uint64_t key = BitPack8(packed);
  
  switch (shape) {
  case Shape8::Pectinate: return lookup8(key, PEC_LOOKUP8);
  case Shape8::Mix:       return lookup8(key, MIX_LOOKUP8);
  case Shape8::Mid:       return lookup8(key, MID_LOOKUP8);
  case Shape8::Balanced:  return lookup8(key, BAL_LOOKUP8);
  }
  return -1;
}

template <size_t N>
int lookup9(uint64_t key, const std::array<SPRScore64, N>& table) {
  auto it = std::lower_bound(
    table.begin(), table.end(), key,
    [](const SPRScore64& a, uint64_t k) {
      return a.key < k;
    }
  );
  
  return (it != table.end() && it->key == key) ? it->score : -1;
}

int lookup9(const SplitSet9& sp1, const SplitSet9& sp2) {
  Shape9 shape = detect_shape9(sp1);
  
  CanonicalInfo9 canon =
    (shape == Shape9::s0) ? canonical9_0(sp1) :
    (shape == Shape9::s1) ? canonical9_1(sp1) :
    (shape == Shape9::s2) ? canonical9_2(sp1) :
    (shape == Shape9::s3) ? canonical9_3(sp1) :
    (shape == Shape9::s4) ? canonical9_4(sp1) :
                            canonical9_5(sp1);
  
  std::array<int,6> packed{};
  for (int i = 0; i < 6; ++i) {
    Split9 s = permute_split9(sp2[i], canon.perm);
    s = polarize9(sp2[i]);
    packed[i] = s;
  }
  
  std::sort(packed.begin(), packed.end());
  uint64_t key = BitPack9(packed);
  
  switch (shape) {
  case Shape9::s0: return lookup9(key, LOOKUP9_0);
  case Shape9::s1: return lookup9(key, LOOKUP9_1);
  case Shape9::s2: return lookup9(key, LOOKUP9_2);
  case Shape9::s3: return lookup9(key, LOOKUP9_3);
  case Shape9::s4: return lookup9(key, LOOKUP9_4);
  case Shape9::s5: return lookup9(key, LOOKUP9_5);
  }
  
  return -1;
}

inline SplitSet6 read_splits6(const Rcpp::RawVector& r) {
  if (r.size() != 3) {
    Rcpp::stop("Expected length-3 raw vector");
  }
  
  SplitSet6 sp{};
  for (int i = 0; i < 3; ++i) {
    sp[i] = static_cast<uint8_t>(r[i]);
  }
  
  return sp;
}

inline SplitSet7 read_splits7(const Rcpp::RawVector& r) {
  if (r.size() != 4) {
    Rcpp::stop("Expected a length-4 raw vector of splits");
  }
  
  SplitSet7 sp{};
  for (int i = 0; i < 4; ++i) {
    sp[i] = static_cast<uint8_t>(r[i]);
  }
  
  return sp;
}

inline SplitSet8 read_splits8(const Rcpp::RawVector& r) {
  if (r.size() != 5)
    Rcpp::stop("Expected length-5 raw vector");
  SplitSet8 sp{};
  for (int i = 0; i < 5; ++i)
    sp[i] = static_cast<uint8_t>(r[i]);
  return sp;
}

inline SplitSet9 read_splits9(const Rcpp::RawVector& r) {
  if (r.size() != 6)
    Rcpp::stop("Expected length-6 raw vector");
  
  SplitSet9 sp{};
  for (int i = 0; i < 6; ++i)
    sp[i] = static_cast<uint16_t>(r[i]);
  
  return sp;
}

// [[Rcpp::export]]
int spr_table_6(const Rcpp::RawVector& sp1, const Rcpp::RawVector& sp2) {
  return lookup6(read_splits6(sp1), read_splits6(sp2));
}

// [[Rcpp::export]]
int spr_table_7(const Rcpp::RawVector& sp1, const Rcpp::RawVector& sp2) {
  SplitSet7 s1 = read_splits7(sp1);
  for (auto& s : s1) s = smaller_split7(s);
  
  return lookup7(s1, read_splits7(sp2));
}

// [[Rcpp::export]]
int spr_table_8(const Rcpp::RawVector& sp1, const Rcpp::RawVector& sp2) {
  SplitSet8 s1 = read_splits8(sp1);
  for (auto& s : s1) s = smaller_split8(s);
 
  return lookup8(s1, read_splits8(sp2));
}

// [[Rcpp::export]]
int spr_table_9(const Rcpp::RawVector& sp1, const Rcpp::RawVector& sp2) {
  SplitSet9 s1 = read_splits9(sp1);
  for (auto& s : s1) s = smaller_split9(s);
  
  return lookup9(s1, read_splits9(sp2));
}
