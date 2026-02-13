#pragma once

using Split6 = uint8_t;                 // lower 6 bits used
using SplitSet6 = std::array<Split6, 3>;
constexpr Split6 MASK6 = 0x3F;          // 00111111

using Split7 = uint8_t;                 // 7 bits used
using SplitSet7 = std::array<Split7, 4>;
using Perm7 = std::array<uint8_t, 7>;
constexpr Split7 MASK7 = 0x7F;          // 01111111

using Split8 = uint8_t;                 // 8 bits used
using SplitSet8 = std::array<Split8, 5>;
using Perm8 = std::array<uint8_t, 8>;
constexpr Split8 MASK8 = 0xFF;

using Split9 = uint16_t;                // 9 bits used
using SplitSet9 = std::array<Split9, 6>;
using Perm9 = std::array<uint16_t, 9>;
constexpr Split9 MASK9 = 0x01FF;

enum class Shape7 { Pectinate, Balanced };
enum class Shape8 { Pectinate, Mix, Mid, Balanced};
enum class Shape9 { s0, s1, s2, s3, s4, s5 };

struct CanonicalInfo7 { Shape7 shape; Perm7 perm; };
struct CanonicalInfo8 { Shape8 shape; Perm8 perm; };
struct CanonicalInfo9 { Shape9 shape; Perm9 perm; };


inline int popcount6(Split6 x) {
  return __builtin_popcount(x & MASK6);
}

inline int popcount7(Split7 x) {
  return __builtin_popcount(x & MASK7);
}

inline int popcount8(Split8 x) {
  return __builtin_popcount(x & MASK8);
}

inline int popcount9(Split9 x) {
  return __builtin_popcount(x & MASK9);
}

inline int tips_in_smallest7(uint8_t x) {
  const int count = popcount7(x);
  return count < 4 ? count : 7 - count;
}

inline int tips_in_smallest8(uint8_t x) {
  const int count = popcount8(x);
  return (count <= 4) ? count : 8 - count;
}

inline int tips_in_smallest9(Split9 s) {
  const int c = popcount9(s);
  return (c <= 4) ? c : 9 - c;
}

inline Split7 xor_split7(Split7 a, Split7 b) {
  return (a ^ b) & MASK7;
}

inline Split8 xor_split8(Split8 a, Split8 b) {
  return (a ^ b) & MASK8;
}

inline Split9 xor_split9(Split9 a, Split9 b) {
  return (a ^ b) & MASK9;
}

inline Split6 smaller_split6(Split6 s) {
  if (popcount6(s) > 3) s ^= MASK6;
  return s;
}

inline Split7 smaller_split7(Split7 s) {
  if (popcount7(s) > 3) s ^= MASK7;
  return s;
}

inline Split8 smaller_split8(Split8 s) {
  if (popcount8(s) > 4) s ^= MASK8;
  return s;
}

inline Split9 smaller_split9(Split9 s) {
  if (popcount9(s) > 4) s ^= MASK9;
  return s;
}

inline Split6 overlapper6(Split6 a, Split6 b) {
  Split6 x = (a ^ b) & MASK6;
  return (popcount6(x) == 1) ? x : (x ^ MASK6);
}
