#ifndef _TREEDIST_INT_H
#define _TREEDIST_INT_H

#include <cstdint>

typedef int_fast16_t int16;
typedef int_fast32_t int32;
typedef uint_fast16_t uint16;
typedef std::vector<typeof(NA_INTEGER)> grf_match;

#define NA_INT16 int16(-0x7fff)
#define UINT_16_MAX uint16(UINT_FAST16_MAX)
#endif