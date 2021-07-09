#ifndef _TREEDIST_INT_H
#define _TREEDIST_INT_H

#include <stdint.h>

typedef int_fast16_t int16;
typedef int_fast32_t int32;
typedef uint_fast16_t uint16;
typedef std::vector<typeof(NA_INTEGER)> grf_match;

const int16 NA_INT16 = -0x7fff;
const int16 INT_16_MAX = int16(INT_FAST16_MAX);
const uint16 UINT_16_MAX = uint16(UINT_FAST16_MAX);
#endif