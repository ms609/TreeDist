#include <Rcpp.h>
#include <algorithm>
#include <array>
#include <cstdint>
#include "spr/lookup_table_7.h"
#include "TreeTools/assert.h"

using Split = uint8_t;          // 7 bits used
using SplitSet = std::array<Split, 4>;
using Perm = std::array<uint8_t, 7>;

inline int popcount7(uint8_t x) {
  return __builtin_popcount(x & 0x7F);
}

inline int tips_in_smallest7(uint8_t x) {
  const int count = __builtin_popcount(x & 0x7F);
  return count < 4 ? count : 7 - count;
}

inline Split xor_split(Split a, Split b) {
  return (a ^ b) & 0x7F;
}

inline Split smaller_split(Split s) {
  if (popcount7(s) > 3) {
    s ^= 0x7F;
  }
  return s;
}

enum class Shape { Pectinate, Balanced };

struct CanonicalInfo {
  Shape shape;
  Perm perm;
};

Shape detect_shape(const SplitSet& sp) {
  int trio_count = 0;
  for (auto s : sp) {
    int t = popcount7(s);
    if (t == 3) ++trio_count;
  }
  return (trio_count == 2) ? Shape::Pectinate : Shape::Balanced;
}

CanonicalInfo canonical_pectinate(const SplitSet& sp) {
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
  
  const Split mid = xor_split(sp[trio1], sp[trio2]) ^ 0x7F;
  
  if (tips_in_smallest7(xor_split(sp[pair1], sp[trio1])) != 1) {
    std::swap(pair1, pair2);
  }
  Split trio1Tip = smaller_split(xor_split(sp[trio1], sp[pair1]));
  Split trio2Tip = smaller_split(xor_split(sp[trio2], sp[pair2]));
  Split duo1 = smaller_split(sp[pair1]);
  Split duo2 = smaller_split(sp[pair2]);
  
  Perm perm{};
  int k = 0;
  
  auto emit = [&](Split s) {
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
  
  return { Shape::Pectinate, perm };
}

CanonicalInfo canonical_balanced(const SplitSet& sp) {
  std::array<int,4> tiss; // Tips in smallest split
  for (int i = 0; i < 4; ++i) {
    tiss[i] = popcount7(sp[i]);
  }
  
  int firstTrio = std::max_element(tiss.begin(), tiss.end()) - tiss.begin();
  Split firstSp = sp[firstTrio];
  
  int trioPair = -1;
  Split solo{};
  for (int i = 0; i < 4; ++i) {
    if (i == firstTrio) continue;
    Split s = xor_split(sp[i], firstSp);
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
  
  Split singleton = solo;
  Split trio = smaller_split(sp[trioPair]);
  Split o1 = smaller_split(sp[other1]);
  Split o2 = smaller_split(sp[other2]);
  
  Perm perm{};
  int k = 0;
  
  auto emit = [&](Split s) {
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
  
  return { Shape::Balanced, perm };
}


Split permute_split(Split s, const Perm& p) {
  Split out = 0;
  for (int i = 0; i < 7; ++i) {
    if (s & (1 << p[i])) {
      out |= (1 << i);
    }
  }
  return out;
}

inline Split polarize(Split s) {
  if (s & (1 << 6)) s ^= 0x7F;
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

int lookup_7(const SplitSet& sp1, const SplitSet& sp2) {
  Shape shape = detect_shape(sp1);
  
  CanonicalInfo canon =
    (shape == Shape::Pectinate)
    ? canonical_pectinate(sp1)
      : canonical_balanced(sp1);
  
  std::array<int,4> packed{};
  for (int i = 0; i < 4; ++i) {
    Split s = permute_split(sp2[i], canon.perm);
    s = polarize(s);
    packed[i] = s;
  }
  
  std::sort(packed.begin(), packed.end());
  
  uint32_t key = BitPack7(packed);
  return (shape == Shape::Pectinate)
    ? lookup(key, PEC_LOOKUP7)
      : lookup(key, BAL_LOOKUP7);
}

inline SplitSet read_splits(const Rcpp::RawVector& r) {
  if (r.size() != 4)
    Rcpp::stop("Expected a length-4 raw vector of splits");
  
  SplitSet sp{};
  for (int i = 0; i < 4; ++i) {
    sp[i] = static_cast<uint8_t>(r[i]);
  }
  
  return sp;
}

inline SplitSet smallest_splits(SplitSet sp) {
  for (auto& s : sp) {
    if (popcount7(s) > 3) {
      s ^= 0x7F;
    }
  }
  return sp;
}

// [[Rcpp::export]]
int spr_table_7(const Rcpp::RawVector& sp1, const Rcpp::RawVector& sp2) {
  SplitSet s1_raw = read_splits(sp1);
  SplitSet s1 = smallest_splits(s1_raw);
  
  SplitSet s2 = read_splits(sp2);
  
  return lookup_7(s1, s2);
}
