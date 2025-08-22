#ifndef _TREEDIST_INT_H
#define _TREEDIST_INT_H

#include <cstdint>

using int16 = int_fast16_t;
using int32 = int_fast32_t;
using uint16 = uint_fast16_t;
using grf_match = std::vector<decltype(NA_INTEGER)> ;

constexpr uint16 INT_16_MAX = INT_FAST16_MAX;
constexpr uint16 UINT_16_MAX = UINT_FAST16_MAX;
constexpr int16 NA_INT16 = -0x7fff;

#endif
