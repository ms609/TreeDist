#ifndef TREEDIST_COL_PRICES_H_
#define TREEDIST_COL_PRICES_H_

#include <cstdint>

// Supplies, on demand, prices for the splits of the second tree in a bounded
// mutual clustering calculation: one per split, in LAP cost units, each in
// [0, BIG].  `get` may return null if no prices are available.  Prices
// affect only speed, never results.
struct ColPrices {
  const int_fast64_t *(*get)(void *context);
  void *context;
};

#endif
