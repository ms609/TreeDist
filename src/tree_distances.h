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
const uint_fast32_t powers_of_two[16] = {1, 2, 4, 8, 16, 32, 64, 128, 256, 512,
                                         1024, 2048, 4096, 8192, 16384, 32768};

const cost ROUND_PRECISION = 2048*2048;

/***** Constants requiring initialization *****/

extern uint_fast32_t bitcounts[65536];
extern double lg2_double_factorial[MAX_TIPS + MAX_TIPS - 2],
  lg2_rooted[MAX_TIPS + 1], lg2_unrooted[MAX_TIPS + 1];

/*************** FUNCTIONS  *******************/

extern int16 count_bits (splitbit x);

extern cost lap(int16 dim, cost **assigncost,
                lap_col *rowsol, lap_row *colsol,
                cost *u, cost *v);

extern double lg2_trees_matching_split(int16 a, int16 b),
  ic_element (const double nkK, const int16 nk,
              const int16 nK, const double n),
  one_overlap (const int16 a, const int16 b, const int16 n),
  one_overlap_notb (const int16 a, const int16 n_minus_b, const int16 n),
  spi (const splitbit* a_state, const splitbit* b_state, const int16 n_tips, 
       const int16 in_a, const int16 in_b, 
       const double lg2_unrooted_n, const int16 n_bins);

extern List cpp_robinson_foulds_distance (RawMatrix x, RawMatrix y, 
                                          IntegerVector nTip);

extern rf_match cpp_robinson_foulds_matching (raw_vector x,
                                              raw_vector y, 
                                              const int16 n_cols,
                                              const int16 n_tips);
