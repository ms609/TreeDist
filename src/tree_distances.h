#include <limits>
#include <Rcpp.h>
#include "ints.h"

using namespace Rcpp;


/*************** TYPES      *******************/

typedef uint_fast64_t splitbit;
typedef int_fast64_t cost;
typedef int16 lap_row;
typedef int16 lap_col;

/*************** CONSTANTS  *******************/

const int16 BIN_SIZE = 64, 
            MAX_BINS = 32,
            MAX_TIPS = BIN_SIZE * MAX_BINS,
            MAX_SPLITS = MAX_TIPS; /* -3, but quicker if a power of two? */

const splitbit ALL_ONES = (std::numeric_limits<splitbit>::max)();

/* For a reason I've not estabilshed, shrinking BIG is necessary to avoid 
 * an infinite loop in lap. */
const cost BIG = ((std::numeric_limits<cost>::max)() / MAX_SPLITS);

const splitbit right16bits = 65535U;
const splitbit powers_of_two[64] = {
  0x1, 0x2, 0x4, 0x8,
  0x10, 0x20, 0x40, 0x80,
  0x100, 0x200, 0x400, 0x800,
  0x1000, 0x2000, 0x4000, 0x8000,
  0x10000, 0x20000, 0x40000, 0x80000,
  0x100000, 0x200000, 0x400000, 0x800000,
  0x1000000, 0x2000000, 0x4000000, 0x8000000,
  0x10000000, 0x20000000, 0x40000000, 0x80000000,
  0x100000000, 0x200000000, 0x400000000, 0x800000000,
  0x1000000000, 0x2000000000, 0x4000000000, 0x8000000000,
  0x10000000000, 0x20000000000, 0x40000000000, 0x80000000000,
  0x100000000000, 0x200000000000, 0x400000000000, 0x800000000000,
  0x1000000000000, 0x2000000000000, 0x4000000000000, 0x8000000000000,
  0x10000000000000, 0x20000000000000, 0x40000000000000, 0x80000000000000,
  0x100000000000000, 0x200000000000000, 0x400000000000000, 0x800000000000000,
  0x1000000000000000, 0x2000000000000000, 0x4000000000000000, 0x8000000000000000
  };

const cost ROUND_PRECISION = 2048 * 2048;

/***** Constants requiring initialization *****/

extern uint_fast32_t bitcounts[65536];
extern double 
  lg2[int32(MAX_TIPS - 1) * (MAX_TIPS - 1) + 1],
  lg2_double_factorial[MAX_TIPS + MAX_TIPS - 2],
  lg2_unrooted[MAX_TIPS + 2];
extern double *lg2_rooted;

/*************** FUNCTIONS  *******************/

extern int16 count_bits (splitbit x);

extern cost lap(int16 dim, cost **assigncost,
                lap_col *rowsol, lap_row *colsol,
                cost *u, cost *v);

extern double 
  mmsi_score(const int16 n_same, const int16 n_a_and_b,
             const int16 n_different, const int16 n_a_only),
  ic_element(const int16 nkK, const int16 nk,
             const int16 nK, const int16 n),
  ic_matching(const int16 a, const int16 b, const int16 n),
  one_overlap(const int16 a, const int16 b, const int16 n),
  one_overlap_notb(const int16 a, const int16 n_minus_b, const int16 n),
  spi_overlap(const splitbit* a_state, const splitbit* b_state, const int16 n_tips, 
      const int16 in_a, const int16 in_b, const int16 n_bins);

extern List cpp_robinson_foulds_distance (RawMatrix x, RawMatrix y, 
                                          IntegerVector nTip);
