#ifndef _TREEDIST_INT_H
#define _TREEDIST_INT_H

#include <cstdint>

typedef int_fast16_t int16;
typedef int_fast32_t int32;
typedef uint_fast16_t uint16;
typedef std::vector<decltype(NA_INTEGER)> grf_match;

constexpr uint16 INT_16_MAX = INT_FAST16_MAX;
constexpr uint16 UINT_16_MAX = UINT_FAST16_MAX;
constexpr int16 NA_INT16 = -0x7fff;

#endif
