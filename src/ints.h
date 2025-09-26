#ifndef _TREEDIST_INT_H
#define _TREEDIST_INT_H

#include <cstdint>

using int16 = int_fast16_t;
using int32 = int_fast32_t;
using uint32 = uint_fast32_t;
using grf_match = std::vector<int>;

constexpr int16 NA_INT16 = std::numeric_limits<int16>::min();
constexpr int32 NA_INT32 = std::numeric_limits<int32>::min();

#endif
