#include <limits>
#include <stdint.h>
#include <Rcpp.h>

using namespace Rcpp;


/*************** TYPES      *******************/

typedef uint64_t splitbit;
typedef int64_t cost;
typedef int lap_row;
typedef int lap_col;

/*************** CONSTANTS  *******************/

const int BIN_SIZE = 64, 
  MAX_BINS = 32,
  MAX_TIPS = BIN_SIZE * MAX_BINS,
  MAX_SPLITS = MAX_TIPS; /* -3, but quicker if a power of two? */

const splitbit ALL_ONES = std::numeric_limits<splitbit>::max();

const double BIGL = double (BIG);
/* For a reason I've not estabilshed, shrinking BIG is necessary to avoid 
 * an infinite loop in lap. */
const cost BIG = (std::numeric_limits<cost>::max() / MAX_SPLITS) / 64;

const splitbit right16bits = 65535U;
const uint32_t powers_of_two[16] = {1, 2, 4, 8, 16, 32, 64, 128, 256, 512,
                                    1024, 2048, 4096, 8192, 16384, 32768};

/***** Constants requiring initialization *****/

extern uint32_t bitcounts[65536];
extern double lg2_double_factorial[MAX_TIPS + MAX_TIPS - 2],
  lg2_rooted[MAX_TIPS + 1], lg2_unrooted[MAX_TIPS + 1];

/*************** FUNCTIONS  *******************/

extern int count_bits (splitbit x);

extern cost lap(int dim, cost **assigncost,
                lap_col *rowsol, lap_row *colsol,
                cost *u, cost *v);

extern double lg2_trees_matching_split(int a, int b),
  p_lg2_p_frac (double p), 
  p_lg2_p (double p),
  entropy2 (double p),
  entropy4 (double p1, double p2, double p3, double p4),
  one_overlap (const int a, const int b, const int n),
  one_overlap_notb (const int a, const int n_minus_b, const int n),
  mpi (const splitbit* a_state, const splitbit* b_state, const int n_tips, 
       const int in_a, const int in_b, 
       const double lg2_unrooted_n, const int n_bins);

extern List cpp_robinson_foulds_distance (RawMatrix x, RawMatrix y, 
                                          IntegerVector nTip);
