#pragma once
#include "lookup.h"
#include "TreeTools/assert.h"

CanonicalInfo8 canonical_pectinate8(const SplitSet8& sp) {
  std::array<int,5> tiss;
  for (int i = 0; i < 5; ++i) {
    tiss[i] = tips_in_smallest8(sp[i]);
  }
  
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
    for (int i = 0; i < 8; ++i) {
      if (s & (1 << i)) {
        perm[k++] = i;
      }
    }
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
