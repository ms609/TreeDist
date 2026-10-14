// Stub for analysis builds: the constants TreeDist/types.h needs, no Rcpp.
#pragma once
#include <cstdint>
using splitbit = uint_fast64_t;
#define R_BIN_SIZE int16_t(8)
#define SL_BIN_SIZE int16_t(64)
#define SL_MAX_BINS int16_t(511)
#define SL_MAX_TIPS (SL_BIN_SIZE * SL_MAX_BINS)
#define SL_MAX_SPLITS (SL_MAX_TIPS - 3)
