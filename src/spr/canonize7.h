#pragma once
#include "lookup.h"
#include "TreeTools/assert.h"

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
  
  const Split7 mid = xor_split7(sp[trio1], sp[trio2]) ^ MASK7;
  
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
      solo = s ^ MASK7;
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
