#ifndef _TREEDIST_LOG2_TABLE_H
#define _TREEDIST_LOG2_TABLE_H

#include <array>
#include <cmath> /* for log2() */
#include <cstdint>

constexpr int_fast32_t LOG_MAX = 2048;

// log2(i) for 0 <= i < LOG_MAX (element 0 is 0), shared by every
// translation unit and filled on first use.
inline const double *log2_lookup_table() {
  static const std::array<double, LOG_MAX> table = [] {
    std::array<double, LOG_MAX> t{};
    for (int i = 1; i < LOG_MAX; ++i) {
      t[i] = std::log2(i);
    }
    return t;
  }();
  return table.data();
}

#endif
