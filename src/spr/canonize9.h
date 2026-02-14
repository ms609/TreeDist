#pragma once
#include "lookup.h"
#include "TreeTools/assert.h"

CanonicalInfo9 canonical9_0(const SplitSet9& sp) {
  std::array<int,6> tiss;
  for (int i = 0; i < 6; ++i)
    tiss[i] = tips_in_smallest9(sp[i]);
  
  std::array<int,2> fours{}, trios{}, pairs{};
  int fi = 0, ti = 0, pi = 0;
  
  for (int i = 0; i < 6; ++i)
  {
    if      (tiss[i] == 4) fours[fi++] = i;
    else if (tiss[i] == 3) trios[ti++] = i;
    else if (tiss[i] == 2) pairs[pi++] = i;
  }
  
  ASSERT(fi == 2 && ti == 2 && pi == 2);
  
  Split9 mid1 = xor_split9(sp[fours[0]], sp[trios[0]]);
  mid1 = smaller_split9(mid1);
  
  Split9 mid2;
  
  if (tips_in_smallest9(mid1) == 1)
  {
    mid2 = smaller_split9(xor_split9(sp[fours[1]], sp[trios[1]]));
  }
  else
  {
    std::swap(trios[0], trios[1]);
    mid1 = smaller_split9(xor_split9(sp[fours[0]], sp[trios[0]]));
    mid2 = smaller_split9(xor_split9(sp[fours[1]], sp[trios[1]]));
  }
  
  if (tips_in_smallest9(xor_split9(sp[trios[0]], sp[pairs[0]])) != 1)
    std::swap(pairs[0], pairs[1]);
  
  Split9 centre = smaller_split9(xor_split9(sp[fours[0]], sp[fours[1]]));
  Split9 trio1  = smaller_split9(xor_split9(sp[trios[0]], sp[pairs[0]]));
  Split9 trio2  = smaller_split9(xor_split9(sp[trios[1]], sp[pairs[1]]));
  
  Split9 pair1 = smaller_split9(sp[pairs[0]]);
  Split9 pair2 = smaller_split9(sp[pairs[1]]);
  
  Perm9 perm{};
  int k = 0;
  
  auto emit = [&](Split9 s)
  {
    for (int i = 0; i < 9; ++i)
      if (s & (1 << i)) perm[k++] = i;
  };
  
  emit(centre);
  emit(mid1);
  emit(trio1);
  emit(pair1);
  emit(mid2);
  emit(trio2);
  emit(pair2);
  
  ASSERT(k == 9);
  return { Shape9::s0, perm };
}

CanonicalInfo9 canonical9_1(const SplitSet9& sp) {
  std::array<int,6> tiss;
  for (int i = 0; i < 6; ++i)
    tiss[i] = tips_in_smallest9(sp[i]);
  
  std::array<int,2> fours{};
  int trio = -1, fi = 0;
  
  for (int i = 0; i < 6; ++i)
  {
    if      (tiss[i] == 4) fours[fi++] = i;
    else if (tiss[i] == 3) trio = i;
  }
  
  ASSERT(fi == 2 && trio >= 0);
  
  Split9 trioSp = sp[trio];
  
  int trioPair = -1;
  Split9 soloSp{};
  
  for (int i = 0; i < 6; ++i)
  {
    if (i == trio || i == fours[0] || i == fours[1]) continue;
    
    Split9 s = xor_split9(trioSp, sp[i]);
    if (tips_in_smallest9(s) == 1)
    {
      trioPair = i;
      soloSp = smaller_split9(s);
      break;
    }
  }
  
  ASSERT(trioPair >= 0);
  
  if (tips_in_smallest9(xor_split9(trioSp, sp[fours[0]])) != 1)
    std::swap(fours[0], fours[1]);
  
  Split9 centre = smaller_split9(xor_split9(sp[fours[0]], trioSp));
  Split9 mid2   = smaller_split9(xor_split9(sp[fours[1]], sp[fours[0]]));
  
  Split9 cherry = smaller_split9(sp[trioPair]);
  
  std::array<Split9,2> others{};
  int oi = 0;
  
  for (int i = 0; i < 6; ++i) {
    if (i != trio && i != fours[0] && i != fours[1] && i != trioPair) {
      others[oi++] = smaller_split9(sp[i]);
    }
  }
  
  Perm9 perm{};
  int k = 0;
  
  auto emit = [&](Split9 s)
  {
    for (int i = 0; i < 9; ++i)
      if (s & (1 << i)) perm[k++] = i;
  };
  
  emit(centre);
  emit(soloSp);
  emit(cherry);
  emit(mid2);
  emit(others[0]);
  emit(others[1]);
  
  ASSERT(k == 9);
  return { Shape9::s1, perm };
}
CanonicalInfo9 canonical9_2(const SplitSet9& sp) {
  std::array<int, 6> tiss{};
  for (int i = 0; i < 6; ++i) {
    tiss[i] = tips_in_smallest9(sp[i]);
  }
  
  // 2) fours = first index where size==4
  int four = -1;
  std::array<int, 2> trios{};
  int ti = 0;
  std::array<int, 3> pairs{};
  int pi = 0;
  
  for (int i = 0; i < 6; ++i) {
    if      (tiss[i] == 4 && four < 0) four = i;
    else if (tiss[i] == 3)             trios[ti++] = i;
  }
  ASSERT(four >= 0 && ti == 2);
  
  // Build 'pairs' in ascending split index, excluding fours and trios
  for (int i = 0; i < 6; ++i) {
    if (i == four || i == trios[0] || i == trios[1]) continue;
    ASSERT(tiss[i] == 2);
    pairs[pi++] = i;
  }
  ASSERT(pi == 3);
  
  // 3) tie-break: if xor(trios[2], fours) is pendant, reverse trios
  {
    Split9 r = xor_split9(sp[trios[1]], sp[four]);
    int c = popcount9(r);
    if (c == 1 || c == 8) std::swap(trios[0], trios[1]);
  }
  
  // 4) trioSp = trios[1], find the first pair that makes pendant xor
  int pair1 = -1, pair2 = -1, pair3 = -1;
  Split9 trioSp  = sp[trios[0]];
  Split9 soloSp1 = 0, soloSp2 = 0;
  
  for (int j = 0; j < 3; ++j) {
    int i = pairs[j];
    Split9 s = xor_split9(trioSp, sp[i]);
    if (tips_in_smallest9(s) == 1) { pair1 = i; soloSp1 = s; break; }
  }
  ASSERT(pair1 >= 0);
  
  // 5) trioSp2 = other trio, find first remaining pair that makes pendant xor
  Split9 trioSp2 = sp[trios[1]];
  for (int j = 0; j < 3; ++j) {
    int i = pairs[j];
    if (i == pair1) continue;
    Split9 s = xor_split9(trioSp2, sp[i]);
    if (tips_in_smallest9(s) == 1) { pair2 = i; soloSp2 = s; break; }
  }
  ASSERT(pair2 >= 0);
  
  // 6) remaining pair
  for (int j = 0; j < 3; ++j) {
    int i = pairs[j];
    if (i != pair1 && i != pair2) { pair3 = i; break; }
  }
  ASSERT(pair3 >= 0);
  
  // 7) Build the six blocks in the exact R order, using "AsTips" semantics
  //    AsTips == smaller side; within a block, tips are in ascending index.
  Split9 block_s = smaller_split9(soloSp1);
  Split9 block_c = smaller_split9(sp[pair1]);
  Split9 block_t = smaller_split9(xor_split9(sp[four], trioSp));
  Split9 block_u = smaller_split9(soloSp2);
  Split9 block_q = smaller_split9(sp[pair2]);
  Split9 block_p = smaller_split9(sp[pair3]);
  
  // 8) Emit FIRST OCCURRENCE ONLY across the concatenation (R's behavior)
  bool seen[9] = {false,false,false,false,false,false,false,false,false};
  Perm9 perm{};
  int k = 0;
  
  auto emit_unique = [&](Split9 s) {
    for (int bit = 0; bit < 9; ++bit) {
      if (s & (Split9(1) << bit)) {
        if (!seen[bit]) {
          seen[bit] = true;
          perm[k++] = bit;
        }
      }
    }
  };
  
  emit_unique(block_s); // s
  emit_unique(block_c); // c
  emit_unique(block_t); // t
  emit_unique(block_u); // u
  emit_unique(block_q); // q
  emit_unique(block_p); // p
  
  ASSERT(k == 9);
  // Verify bijection
  for (int i = 0; i < 9; ++i) ASSERT(seen[i]);
  
  return { Shape9::s2, perm };
}
inline int single_tip(Split9 s) {
  ASSERT(s == (s &= MASK9));
  const int c = popcount9(s);
  if (c == 1) return __builtin_ctz(s);
  if (c == 8) return __builtin_ctz((~s) & MASK9);
  Rcpp::Rcout << int(s);
  Rcpp::stop("single_tip(): split is not singleton-sized");
}

inline Split9 polarize9(Split9 s, int tip) {
  s &= MASK9;
  
  if (!(s & (Split9(1) << tip))) {
    s ^= MASK9;
  }
  
  return s;
}

CanonicalInfo9 canonical9_3(const SplitSet9& sp) {
  std::array<int, 6> tiss{};
  for (int i = 0; i < 6; ++i) {
    tiss[i] = tips_in_smallest9(sp[i]);
  }
  
  std::array<int, 2> fours{};
  std::array<int, 1> trios{};
  std::array<int, 3> pairs{};
  
  int fi = 0, ti = 0, pi = 0;
  
  for (int i = 0; i < 6; ++i) {
    if (tiss[i] == 4) {
      fours[fi++] = i;
    } else if (tiss[i] == 3) {
      trios[ti++] = i;
    } else if (tiss[i] == 2) {
      pairs[pi++] = i;
    }
  }
  
  ASSERT(fi == 2 && ti == 1 && pi == 3);
  
  const int trio = trios[0];
  Split9 trioSp = sp[trio];
  
  int trioPair = -1;
  Split9 soloSp{};
  
  for (int i = 0; i < 3; ++i) {
    const int p = pairs[i];
    
    Split9 s = xor_split9(trioSp, sp[p]);
    if (tips_in_smallest9(s) == 1) {
      trioPair = p;
      soloSp = smaller_split9(s);
      break;
    }
  }
  
  ASSERT(trioPair >= 0);
  
  const int soloTip = single_tip(soloSp);
  
  Split9 quadSp = polarize9(sp[fours[0]], soloTip) &
    polarize9(sp[fours[1]], soloTip);
  
  int midPair = -1;
  
  for (int i = 0; i < 3; ++i) {
    const int p = pairs[i];
    if (p == trioPair) {
      continue;
    }
    
    Split9 mid = quadSp & polarize9(sp[p], soloTip);
    if (tips_in_smallest9(mid) == 3) {
      midPair = p;
      break;
    }
  }
  
  ASSERT(midPair >= 0);
  
  std::array<int, 2> remaining{};
  int ri = 0;
  
  for (int i = 0; i < 3; ++i) {
    const int p = pairs[i];
    if (p != trioPair && p != midPair) {
      remaining[ri++] = p;
    }
  }
  
  ASSERT(ri == 2);
  
  Perm9 perm{};
  int k = 0;
  
  auto emit = [&](Split9 s) {
    s = smaller_split9(s);
    for (int i = 0; i < 9; ++i) {
      if (s & (Split9(1) << i)) {
        perm[k++] = i;
      }
    }
  };
  
  emit(soloSp);
  emit(sp[trioPair]);
  emit(sp[midPair]);
  emit(sp[remaining[0]]);
  emit(sp[remaining[1]]);
  
  ASSERT(k == 9);
  return { Shape9::s3, perm };
}

CanonicalInfo9 canonical9_4(const SplitSet9& sp) {
  std::array<int, 6> tiss{};
  for (int i = 0; i < 6; ++i) {
    tiss[i] = tips_in_smallest9(sp[i]);
  }
  
  std::array<int, 3> trios{};
  std::array<int, 3> pairs{};
  
  int ti = 0, pi = 0;
  
  for (int i = 0; i < 6; ++i) {
    if (tiss[i] == 3) {
      trios[ti++] = i;
    } else if (tiss[i] == 2) {
      pairs[pi++] = i;
    }
  }
  
  ASSERT(ti == 3 && pi == 3);
  
  std::array<Split9, 3> solo{};
  std::array<int, 3> pairIdx{};
  
  std::array<bool, 3> usedPair{ false, false, false };
  
  for (int t = 0; t < 3; ++t) {
    Split9 trioSp = sp[trios[t]];
    
    for (int p = 0; p < 3; ++p) {
      if (usedPair[p]) {
        continue;
      }
      
      Split9 s = xor_split9(trioSp, sp[pairs[p]]);
      if (tips_in_smallest9(s) == 1) {
        solo[t] = smaller_split9(s);
        pairIdx[t] = pairs[p];
        usedPair[p] = true;
        break;
      }
    }
  }
  
  Perm9 perm{};
  int k = 0;
  
  auto emit = [&](Split9 s) {
    s = smaller_split9(s);
    for (int i = 0; i < 9; ++i) {
      if (s & (Split9(1) << i)) {
        perm[k++] = i;
      }
    }
  };
  
  for (int i = 0; i < 3; ++i) {
    emit(solo[i]);
    emit(sp[pairIdx[i]]);
  }
  
  ASSERT(k == 9);
  return { Shape9::s4, perm };
}

CanonicalInfo9 canonical9_5(const SplitSet9& sp) {
  std::array<int, 6> tiss{};
  for (int i = 0; i < 6; ++i) {
    tiss[i] = tips_in_smallest9(sp[i]);
  }
  
  std::array<int, 2> fours{};
  std::array<int, 4> pairs{};
  
  int fi = 0, pi = 0;
  
  for (int i = 0; i < 6; ++i) {
    if (tiss[i] == 4) {
      fours[fi++] = i;
    } else if (tiss[i] == 2) {
      pairs[pi++] = i;
    }
  }
  
  ASSERT(fi == 2 && pi == 4);
  
  Split9 soloSp = smaller_split9(xor_split9(sp[fours[0]], sp[fours[1]]));
  const int soloTip = single_tip(soloSp);
  
  Split9 four1 = polarize9(sp[fours[0]], soloTip);
  
  std::array<int, 2> cluster1{};
  std::array<int, 2> cluster2{};
  int c1 = 0, c2 = 0;
  
  for (int i = 0; i < 4; ++i) {
    const int p = pairs[i];
    
    Split9 s = polarize9(sp[p], soloTip);
    
    if (four1 & s) {
      cluster1[c1++] = p;
    } else {
      cluster2[c2++] = p;
    }
  }
  
  ASSERT(c1 == 2 && c2 == 2);
  
  Perm9 perm{};
  int k = 0;
  
  auto emit = [&](Split9 s) {
    s = smaller_split9(s);
    for (int i = 0; i < 9; ++i) {
      if (s & (Split9(1) << i)) {
        perm[k++] = i;
      }
    }
  };
  
  emit(soloSp);
  emit(sp[cluster1[0]]);
  emit(sp[cluster1[1]]);
  emit(sp[cluster2[0]]);
  emit(sp[cluster2[1]]);
  
  ASSERT(k == 9);
  return { Shape9::s5, perm };
}
