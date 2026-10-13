#ifndef TREEDIST_POPCOUNT_H_
#define TREEDIST_POPCOUNT_H_

// Population count that is safe on x86-64 CPUs without POPCNT
// (ms609/TreeTools#287), at full speed where POPCNT exists.
//
// TreeTools versions that define TREETOOLS_POPCOUNT_DISPATCH provide these
// definitions; with older TreeTools, identical local copies are used, so
// TreeDist builds against either.
//
// - count_bits(x): drop-in; one well-predicted branch on a load-time CPUID
//   flag per call.
// - dispatch_popcount(f): calls f(PopcountHW{}) or f(PopcountSW{}) once, so
//   a kernel templated on the policy and calling Pop::count(x) runs exactly
//   as fast as unconditional POPCNT.
// - PopcountDefault: policy that calls count_bits(); the default for helpers
//   that are also called from outside a dispatched kernel.
//
// Note: with older TreeTools, TreeTools' own SplitList constructor still
// executes POPCNT unconditionally; only a TreeTools update removes that.

#include <cstdint>
#include <TreeTools/SplitList.h>

namespace TreeDist {

#ifdef TREETOOLS_POPCOUNT_DISPATCH

  using TreeTools::count_bits;
  using TreeTools::dispatch_popcount;
  using TreeTools::has_popcnt;
  using TreeTools::PopcountHW;
  using TreeTools::PopcountSW;

#else

  inline int32_t popcount_sw(splitbit x) noexcept {
#if defined(__GNUC__) || defined(__clang__)
    return static_cast<int32_t>(__builtin_popcountll(x));
#else
    int32_t count = 0;
    while (x != 0) {
      x &= (x - 1);
      ++count;
    }
    return count;
#endif
  }

#if (defined(__GNUC__) || defined(__clang__)) && defined(__x86_64__) && \
    !defined(__POPCNT__)
  inline bool cpu_has_popcnt() noexcept {
    __builtin_cpu_init();
    return __builtin_cpu_supports("popcnt");
  }

  // Evaluated once, when the shared object is loaded. Anything that runs
  // before dynamic initialization reads false and takes the portable path.
  inline const bool has_popcnt = cpu_has_popcnt();

  // Executes POPCNT unconditionally: call only where has_popcnt is true.
  // volatile stops the compiler hoisting the instruction above that check.
  inline int32_t popcount_hw(splitbit x) noexcept {
    uint64_t result;
    __asm__ __volatile__ ("popcnt %1, %0" : "=r" (result) : "r" (x));
    return static_cast<int32_t>(result);
  }
#elif defined(_MSC_VER) && defined(_M_X64)
} // namespace TreeDist
#include <intrin.h>
namespace TreeDist {
  inline constexpr bool has_popcnt = true;
  inline int32_t popcount_hw(splitbit x) noexcept {
    return static_cast<int32_t>(__popcnt64(x));
  }
#else
  // Non-x86 (the builtin maps to the native instruction), or x86-64 built
  // with -mpopcnt / -march=x86-64-v2 or later (it compiles to POPCNT).
  inline constexpr bool has_popcnt = true;
  inline int32_t popcount_hw(splitbit x) noexcept { return popcount_sw(x); }
#endif

  inline int32_t count_bits(splitbit x) noexcept {
#if defined(__GNUC__) || defined(__clang__)
    if (__builtin_expect(has_popcnt, 1)) return popcount_hw(x);
#else
    if (has_popcnt) return popcount_hw(x);
#endif
    return popcount_sw(x);
  }

  struct PopcountHW {
    static int32_t count(splitbit x) noexcept { return popcount_hw(x); }
  };
  struct PopcountSW {
    static int32_t count(splitbit x) noexcept { return popcount_sw(x); }
  };

  template <typename F>
  inline decltype(auto) dispatch_popcount(F&& f) {
    if (has_popcnt) return f(PopcountHW{});
    return f(PopcountSW{});
  }

#endif // TREETOOLS_POPCOUNT_DISPATCH

  struct PopcountDefault {
    static auto count(splitbit x) noexcept { return count_bits(x); }
  };

} // namespace TreeDist

#endif // TREEDIST_POPCOUNT_H_
