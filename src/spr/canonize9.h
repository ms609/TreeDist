#pragma once
#include "lookup.h"
#include "TreeTools/assert.h"

CanonicalInfo9 canonical9_0(const SplitSet9& sp)
{
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
CanonicalInfo9 canonical9_1(const SplitSet9& sp)
{
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

CanonicalInfo9 canonical9_2(const SplitSet9& sp)
{
  std::array<int,6> tiss;
  for (int i = 0; i < 6; ++i) {
    tiss[i] = tips_in_smallest9(sp[i]);
  }
  
  int four = -1;
  std::array<int,2> trios{};
  int ti = 0;
  
  for (int i = 0; i < 6; ++i) {
    if      (tiss[i] == 4) four = i;
    else if (tiss[i] == 3) trios[ti++] = i;
  }
  
  ASSERT(four >= 0 && ti == 2);
  
  Split9 trioSp1 = sp[trios[0]];
  Split9 trioSp2 = sp[trios[1]];
  
  Split9 solo1{}, solo2{};
  int pair1 = -1, pair2 = -1;
  
  for (int i = 0; i < 6; ++i) {
    if (i == four || i == trios[0] || i == trios[1]) continue;
    
    Split9 s1 = xor_split9(trioSp1, sp[i]);
    if (tips_in_smallest9(s1) == 1)
    {
      pair1 = i;
      solo1 = smaller_split9(s1);
      break;
    }
  }
  
  for (int i = 0; i < 6; ++i)
  {
    if (i == four || i == trios[0] || i == trios[1] || i == pair1) continue;
    
    Split9 s2 = xor_split9(trioSp2, sp[i]);
    if (tips_in_smallest9(s2) == 1)
    {
      pair2 = i;
      solo2 = smaller_split9(s2);
      break;
    }
  }
  
  ASSERT(pair1 >= 0 && pair2 >= 0);
  
  Split9 t = smaller_split9(xor_split9(sp[four], trioSp1));
  
  Perm9 perm{};
  int k = 0;
  
  auto emit = [&](Split9 s)
  {
    for (int i = 0; i < 9; ++i)
      if (s & (1 << i)) perm[k++] = i;
  };
  
  emit(solo1);
  emit(smaller_split9(sp[pair1]));
  emit(t);
  emit(solo2);
  emit(smaller_split9(sp[pair2]));
  
  ASSERT(k == 9);
  return { Shape9::s2, perm };
}

CanonicalInfo9 canonical9_3(const SplitSet9& sp)
{
  std::array<int, 6> tiss{};
  for (int i = 0; i < 6; ++i) {
    int k = popcount9(sp[i]);
    tiss[i] = std::min(k, 9 - k);
  }
  
  int trio = -1, four = -1;
  for (int i = 0; i < 6; ++i) {
    if (tiss[i] == 3 && trio < 0) trio = i;
    if (tiss[i] == 4 && four < 0) four = i;
  }
  
  Split9 trioSp = sp[trio];
  
  int trioPair = -1;
  Split9 soloSp{};
  
  for (int i = 0; i < 6; ++i) {
    if (i == trio || tiss[i] != 2) continue;
    
    Split9 x = xor_split9(trioSp, sp[i]);
    int k = popcount9(x);
    
    if (k == 1 || k == 8) {
      trioPair = i;
      soloSp = x;
      break;
    }
  }
  
  int soloTip = single_tip(soloSp);
  
  Split9 fourSp = polarize9(sp[four], soloTip);
  
  int midPair = -1;
  for (int i = 0; i < 6; ++i) {
    if (i == trio || i == four || i == trioPair || tiss[i] != 2)
      continue;
    
    Split9 mid = fourSp & polarize9(sp[i], soloTip);
    if (popcount9(mid) == 3) {
      midPair = i;
      break;
    }
  }
  
  std::vector<int> remaining;
  for (int i = 0; i < 6; ++i)
    if (tiss[i] == 2 && i != trioPair && i != midPair)
      remaining.push_back(i);
    
    CanonicalInfo9 out{};
    
    out.order = {
      soloTip,
      single_tip(sp[trioPair]),
      single_tip(sp[midPair]),
      single_tip(sp[remaining[0]]),
      single_tip(sp[remaining[1]]),
      -1, -1, -1, -1
    };
    
    return out;
}
CanonicalInfo9 canonical9_4(const SplitSet9& sp)
{
  std::array<int, 6> tiss{};
  for (int i = 0; i < 6; ++i) {
    int k = popcount9(sp[i]);
    tiss[i] = std::min(k, 9 - k);
  }
  
  std::vector<int> trios, pairs;
  for (int i = 0; i < 6; ++i) {
    if (tiss[i] == 3) trios.push_back(i);
    else if (tiss[i] == 2) pairs.push_back(i);
  }
  
  int pair1 = -1, pair2 = -1;
  
  Split9 solo1{}, solo2{}, solo3{};
  
  for (int p : pairs) {
    Split9 x = xor_split9(sp[trios[0]], sp[p]);
    int k = popcount9(x);
    if (k == 1 || k == 8) {
      pair1 = p;
      solo1 = x;
      break;
    }
  }
  
  for (int p : pairs) {
    if (p == pair1) continue;
    Split9 x = xor_split9(sp[trios[1]], sp[p]);
    int k = popcount9(x);
    if (k == 1 || k == 8) {
      pair2 = p;
      solo2 = x;
      break;
    }
  }
  
  int pair3 = -1;
  for (int p : pairs)
    if (p != pair1 && p != pair2)
      pair3 = p;
    
    solo3 = xor_split9(sp[trios[2]], sp[pair3]);
    
    CanonicalInfo9 out{};
    
    out.order = {
      single_tip(solo1),
      single_tip(sp[pair1]),
      single_tip(solo2),
      single_tip(sp[pair2]),
      single_tip(solo3),
      single_tip(sp[pair3]),
      -1, -1, -1
    };
    
    return out;
}

CanonicalInfo9 canonical9_5(const SplitSet9& sp)
{
  std::array<int, 6> tiss{};
  for (int i = 0; i < 6; ++i) {
    int k = popcount9(sp[i]);
    tiss[i] = std::min(k, 9 - k);
  }
  
  std::vector<int> fours, pairs;
  for (int i = 0; i < 6; ++i) {
    if (tiss[i] == 4) fours.push_back(i);
    else if (tiss[i] == 2) pairs.push_back(i);
  }
  
  Split9 soloSp = xor_split9(sp[fours[0]], sp[fours[1]]);
  int soloTip = single_tip(soloSp);
  
  Split9 four1 = polarize9(sp[fours[0]], soloTip);
  
  std::vector<int> group1, group2;
  
  for (int p : pairs) {
    Split9 x = four1 & polarize9(sp[p], soloTip);
    if (popcount9(x) > 0) group1.push_back(p);
    else group2.push_back(p);
  }
  
  CanonicalInfo9 out{};
  
  out.order = {
    soloTip,
    single_tip(sp[group1[0]]),
    single_tip(sp[group1[1]]),
    single_tip(sp[group2[0]]),
    single_tip(sp[group2[1]]),
    -1, -1, -1, -1
  };
  
  return out;
}
