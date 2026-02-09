#include <Rcpp.h>
#include <algorithm>
#include <array>
#include <cstdint>
#include "spr/lookup_table_7.h"
#include "spr/lookup_table_8.h"
#include "TreeTools/assert.h"

using Split7 = uint8_t;          // 7 bits used
using SplitSet7 = std::array<Split7, 4>;
using Perm7 = std::array<uint8_t, 7>;

using Split8 = uint8_t;                 // 8 bits used
using SplitSet8 = std::array<Split8, 5>;
using Perm8 = std::array<uint8_t, 8>;

using Split9 = uint16_t;                 // 9 bits used
using SplitSet9 = std::array<Split9, 6>;
using Perm9 = std::array<uint16_t, 9>;

inline int popcount7(uint8_t x) {
  return __builtin_popcount(x & 0x7F);
}

inline int popcount8(uint8_t x) {
  return __builtin_popcount(x);
}

inline int tips_in_smallest7(uint8_t x) {
  const int count = __builtin_popcount(x & 0x7F);
  return count < 4 ? count : 7 - count;
}

inline int tips_in_smallest8(uint8_t x) {
  int c = popcount8(x);
  return (c <= 4) ? c : 8 - c;
}

inline Split7 xor_split7(Split7 a, Split7 b) {
  return (a ^ b) & 0x7F;
}

inline Split8 xor_split8(Split8 a, Split8 b) {
  return a ^ b;
}

inline Split7 smaller_split7(Split7 s) {
  if (popcount7(s) > 3) {
    s ^= 0x7F;
  }
  return s;
}

inline Split8 smaller_split8(Split8 s) {
  if (popcount8(s) > 4) s ^= 0xFF;
  return s;
}

enum class Shape7 { Pectinate, Balanced };
enum class Shape8 { Pectinate, Mix, Mid, Balanced};
enum class Shape9 { s0, s1, s2, s3, s4, s5 };


struct CanonicalInfo7 {
  Shape7 shape;
  Perm7 perm;
};

struct CanonicalInfo8 {
  Shape8 shape;
  Perm8 perm;
};

struct CanonicalInfo9 {
  Shape9 shape;
  Perm9 perm;
};

inline Shape7 detect_shape7(const SplitSet7& sp) {
  int trio_count = 0;
  for (auto s : sp) {
    int t = popcount7(s);
    if (t == 3) ++trio_count;
  }
  return (trio_count == 2) ? Shape7::Pectinate : Shape7::Balanced;
}

inline Shape8 detect_shape8(const SplitSet8& sp) {
  int n4 = 0, n3 = 0, n2 = 0;
  for (auto s : sp) {
    int t = tips_in_smallest8(s);
    if (t == 4) ++n4;
    else if (t == 3) ++n3;
    else if (t == 2) ++n2;
  }
  
  if (n4 == 1 && n3 == 2) return Shape8::Pectinate;
  if (n4 == 1 && n3 == 1) return Shape8::Mix;
  if (n3 == 2)           return Shape8::Mid;
  return Shape8::Balanced;
}

CanonicalInfo7 canonical_pectinate(const SplitSet7& sp) {
  std::array<int,4> tiss; // Tips in split (smallest)
  for (int i = 0; i < 4; ++i) {
    tiss[i] = popcount7(sp[i]);
  }
  
  int trio1 = -1, trio2 = -1;
  int p = 0;
  for (int i = 0; i < 4; ++i) {
    if (tiss[i] == 3) {
      (p++ == 0 ? trio1 : trio2) = i;
    }
  }
  int pair1 = -1, pair2 = -1;
  for (int i = 0; i < 4; ++i) {
    if (i != trio1 && i != trio2) {
      (pair1 == -1 ? pair1 : pair2) = i;
    }
  }
  
  const Split7 mid = xor_split7(sp[trio1], sp[trio2]) ^ 0x7F;
  
  if (tips_in_smallest7(xor_split7(sp[pair1], sp[trio1])) != 1) {
    std::swap(pair1, pair2);
  }
  Split7 trio1Tip = smaller_split7(xor_split7(sp[trio1], sp[pair1]));
  Split7 trio2Tip = smaller_split7(xor_split7(sp[trio2], sp[pair2]));
  Split7 duo1 = smaller_split7(sp[pair1]);
  Split7 duo2 = smaller_split7(sp[pair2]);
  
  Perm7 perm{};
  int k = 0;
  
  auto emit = [&](Split7 s) {
    for (int i = 0; i < 7; ++i) {
      if (s & (1 << i)) {
        perm[k++] = i;
      }
    }
  };
  
  emit(mid);
  emit(trio1Tip);
  emit(duo1);
  emit(trio2Tip);
  emit(duo2);
  
  return { Shape7::Pectinate, perm };
}

CanonicalInfo7 canonical_balanced(const SplitSet7& sp) {
  std::array<int,4> tiss; // Tips in smallest split
  for (int i = 0; i < 4; ++i) {
    tiss[i] = popcount7(sp[i]);
  }
  
  int firstTrio = std::max_element(tiss.begin(), tiss.end()) - tiss.begin();
  Split7 firstSp = sp[firstTrio];
  
  int trioPair = -1;
  Split7 solo{};
  for (int i = 0; i < 4; ++i) {
    if (i == firstTrio) continue;
    Split7 s = xor_split7(sp[i], firstSp);
    const int s_count = popcount7(s);
    if (s_count == 1) {
      trioPair = i;
      solo = s;
      break;
    } else if (s_count == 6) {
      trioPair = i;
      solo = s ^ 0x7F;
      break;
    }
  }
  ASSERT(trioPair > -1);
  
  int other1 = -1, other2 = -1;
  for (int i = 0; i < 4; ++i) {
    if (i != firstTrio && i != trioPair)
      (other1 == -1 ? other1 : other2) = i;
  }
  ASSERT(other1 > -1);
  ASSERT(other2 > -1);
  
  Split7 singleton = solo;
  Split7 trio = smaller_split7(sp[trioPair]);
  Split7 o1 = smaller_split7(sp[other1]);
  Split7 o2 = smaller_split7(sp[other2]);
  
  Perm7 perm{};
  int k = 0;
  
  auto emit = [&](Split7 s) {
    for (int i = 0; i < 7; ++i) {
      if (s & (1 << i)) {
        perm[k++] = i;
      }
    }
  };
  
  emit(singleton);
  emit(trio);
  emit(o1);
  emit(o2);
  
  return { Shape7::Balanced, perm };
}

CanonicalInfo8 canonical_pectinate8(const SplitSet8& sp) {
  std::array<int,5> tiss;
  for (int i = 0; i < 5; ++i)
    tiss[i] = tips_in_smallest8(sp[i]);
  
  int quad = -1;
  std::array<int,2> trios{};
  std::array<int,2> pairs{};
  int ti = 0, pi = 0;
  
  for (int i = 0; i < 5; ++i) {
    if (tiss[i] == 4) quad = i;
    else if (tiss[i] == 3) trios[ti++] = i;
    else if (tiss[i] == 2) pairs[pi++] = i;
  }
  
  ASSERT(quad >= 0 && ti == 2 && pi == 2);
  
  Split8 mid1 = xor_split8(sp[quad], sp[trios[0]]);
  Split8 mid2 = xor_split8(sp[quad], sp[trios[1]]);
  
  mid1 = smaller_split8(mid1);
  mid2 = smaller_split8(mid2);
  
  // Align pairs to trios
  if (tips_in_smallest8(xor_split8(sp[trios[0]], sp[pairs[0]])) != 1)
    std::swap(pairs[0], pairs[1]);
  
  Split8 trio1Tip = smaller_split8(xor_split8(sp[trios[0]], sp[pairs[0]]));
  Split8 trio2Tip = smaller_split8(xor_split8(sp[trios[1]], sp[pairs[1]]));
  
  Split8 duo1 = smaller_split8(sp[pairs[0]]);
  Split8 duo2 = smaller_split8(sp[pairs[1]]);
  
  Perm8 perm{};
  int k = 0;
  
  auto emit = [&](Split8 s) {
    for (int i = 0; i < 8; ++i)
      if (s & (1 << i)) perm[k++] = i;
  };
  
  emit(mid1);
  emit(trio1Tip);
  emit(duo1);
  emit(mid2);
  emit(trio2Tip);
  emit(duo2);
  
  ASSERT(k == 8);
  return { Shape8::Pectinate, perm };
}
CanonicalInfo8 canonical_mix8(const SplitSet8& sp) {
  std::array<int,5> tiss;
  for (int i = 0; i < 5; ++i)
    tiss[i] = tips_in_smallest8(sp[i]);
  
  int quad = -1, trio = -1;
  for (int i = 0; i < 5; ++i) {
    if (tiss[i] == 4) quad = i;
    else if (tiss[i] == 3) trio = i;
  }
  ASSERT(quad >= 0 && trio >= 0);
  
  Split8 trioSp = sp[trio];
  Split8 mid = xor_split8(sp[quad], trioSp);
  mid = smaller_split8(mid);
  
  int trioPair = -1;
  for (int i = 0; i < 5; ++i) {
    if (i == quad || i == trio) continue;
    Split8 solo = xor_split8(trioSp, sp[i]);
    if (tips_in_smallest8(solo) == 1) {
      trioPair = i;
      break;
    }
  }
  ASSERT(trioPair >= 0);
  
  Split8 trioTip = smaller_split8(xor_split8(trioSp, sp[trioPair]));
  Split8 trioPairTip = smaller_split8(sp[trioPair]);
  
  std::array<Split8,2> others{};
  int oi = 0;
  for (int i = 0; i < 5; ++i) {
    if (i != quad && i != trio && i != trioPair) {
      others[oi++] = smaller_split8(sp[i]);
    }
  }
  
  Perm8 perm{};
  int k = 0;
  
  auto emit = [&](Split8 s) {
    for (int i = 0; i < 8; ++i) {
      if (s & (1 << i)) perm[k++] = i;
    }
  };
  
  emit(mid);
  emit(trioTip);
  emit(trioPairTip);
  emit(others[0]);
  emit(others[1]);
  
  ASSERT(k == 8);
  return { Shape8::Mix, perm };
}
CanonicalInfo8 canonical_mid8(const SplitSet8& sp) {
  std::array<int,5> tiss;
  for (int i = 0; i < 5; ++i)
    tiss[i] = tips_in_smallest8(sp[i]);
  
  std::array<int,2> trios{};
  int ti = 0;
  for (int i = 0; i < 5; ++i) {
    if (tiss[i] == 3) trios[ti++] = i;
  }
    
  ASSERT(ti == 2);
  
  auto find_pair = [&](int trio) {
    for (int i = 0; i < 5; ++i) {
      if (i == trios[0] || i == trios[1]) continue;
      Split8 solo = xor_split8(sp[i], sp[trio]);
      if (tips_in_smallest8(solo) == 1) return i;
    }
    return -1;
  };
  
  int p1 = find_pair(trios[0]);
  int p2 = find_pair(trios[1]);
  ASSERT(p1 >= 0 && p2 >= 0 && p1 != p2);
  
  Split8 solo1 = smaller_split8(xor_split8(sp[p1], sp[trios[0]]));
  Split8 duo1  = smaller_split8(sp[p1]);
  Split8 solo2 = smaller_split8(xor_split8(sp[p2], sp[trios[1]]));
  Split8 duo2  = smaller_split8(sp[p2]);
  
  int rem = 0;
  for (int i = 0; i < 5; ++i) {
    if (i != trios[0] && i != trios[1] && i != p1 && i != p2) {
      rem = i;
    }
  }
  
  Split8 last = smaller_split8(sp[rem]);
  
  Perm8 perm{};
  int k = 0;
  
  auto emit = [&](Split8 s) {
    for (int i = 0; i < 8; ++i) {
      if (s & (1 << i)) perm[k++] = i;
    }
  };
  
  emit(solo1);
  emit(duo1);
  emit(solo2);
  emit(duo2);
  emit(last);
  
  ASSERT(k == 8);
  return { Shape8::Mid, perm };
}

CanonicalInfo8 canonical_balanced8(const SplitSet8& sp) {
  Perm8 perm{};
  int k = 0;
  
  for (auto s : sp) {
    if (tips_in_smallest8(s) == 2) {
      Split8 c = smaller_split8(s);
      for (int i = 0; i < 8; ++i)
        if (c & (1 << i)) perm[k++] = i;
    }
  }
  
  ASSERT(k == 8);
  return { Shape8::Balanced, perm };
}

Split7 permute_split(Split7 s, const Perm7& p) {
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
    if (s & (1 << p[i])) out |= (1 << i);
  }
  return out;
}

inline Split7 polarize7(Split7 s) {
  if (s & (1 << 6)) s ^= 0x7F;
  return s;
}

inline Split8 polarize8(Split8 s) {
  if (s & (1 << 7)) s ^= 0xFF;
  return s;
}

inline uint32_t BitPack7(const std::array<int,4>& v) {
  return ((v[0] - 3)  << 18) |
    ((v[1] - 7)  << 12) |
    ((v[2] - 15) << 6)  |
    ( v[3] - 33);
}

template <size_t N>
int lookup(uint32_t key, const std::array<SPRScore, N>& table) {
  auto it = std::lower_bound(
    table.begin(), table.end(), key,
    [](const SPRScore& a, uint32_t k) { return a.key < k; }
  );
  return (it != table.end() && it->key == key) ? it->score : -1;
}

inline uint64_t BitPack8(const std::array<int,5>& v) {
  return
  ((uint64_t)(v[0]) << 27) |
    ((uint64_t)(v[1]) << 20) |
    ((uint64_t)(v[2]) << 13) |
    ((uint64_t)(v[3]) << 6)  |
    (uint64_t)(v[4]);
}

template <size_t N>
int lookup8(uint64_t key, const std::array<SPRScore64, N>& table) {
  auto it = std::lower_bound(
    table.begin(), table.end(), key,
    [](const SPRScore64& a, uint64_t k) { return a.key < k; }
  );
  return (it != table.end() && it->key == key) ? it->score : -1;
}

int lookup_7(const SplitSet7& sp1, const SplitSet7& sp2) {
  Shape7 shape = detect_shape7(sp1);
  
  CanonicalInfo7 canon =
    (shape == Shape7::Pectinate)
    ? canonical_pectinate(sp1)
      : canonical_balanced(sp1);
  
  std::array<int,4> packed{};
  for (int i = 0; i < 4; ++i) {
    Split7 s = permute_split(sp2[i], canon.perm);
    s = polarize7(s);
    packed[i] = s;
  }
  
  std::sort(packed.begin(), packed.end());
  
  uint32_t key = BitPack7(packed);
  return (shape == Shape7::Pectinate)
    ? lookup(key, PEC_LOOKUP7)
      : lookup(key, BAL_LOOKUP7);
}
int lookup_8(const SplitSet8& sp1, const SplitSet8& sp2) {
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

inline SplitSet7 read_splits(const Rcpp::RawVector& r) {
  if (r.size() != 4)
    Rcpp::stop("Expected a length-4 raw vector of splits");
  
  SplitSet7 sp{};
  for (int i = 0; i < 4; ++i) {
    sp[i] = static_cast<uint8_t>(r[i]);
  }
  
  return sp;
}

inline SplitSet7 smallest_splits(SplitSet7 sp) {
  for (auto& s : sp) {
    if (popcount7(s) > 3) {
      s ^= 0x7F;
    }
  }
  return sp;
}

// [[Rcpp::export]]
int spr_table_7(const Rcpp::RawVector& sp1, const Rcpp::RawVector& sp2) {
  SplitSet7 s1_raw = read_splits(sp1);
  SplitSet7 s1 = smallest_splits(s1_raw);
  
  SplitSet7 s2 = read_splits(sp2);
  
  return lookup_7(s1, s2);
}

inline SplitSet8 read_splits8(const Rcpp::RawVector& r) {
  if (r.size() != 5)
    Rcpp::stop("Expected length-5 raw vector");
  SplitSet8 sp{};
  for (int i = 0; i < 5; ++i)
    sp[i] = static_cast<uint8_t>(r[i]);
  return sp;
}

// [[Rcpp::export]]
int spr_table_8(const Rcpp::RawVector& sp1,
                const Rcpp::RawVector& sp2) {
  SplitSet8 s1 = read_splits8(sp1);
  for (auto& s : s1) s = smaller_split8(s);
  
  SplitSet8 s2 = read_splits8(sp2);
  return lookup_8(s1, s2);
}
