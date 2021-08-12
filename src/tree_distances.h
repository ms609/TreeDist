#ifndef _TREEDIST_TREE_DISTANCES_H
#define _TREEDIST_TREE_DISTANCES_H

#include <TreeTools/SplitList.h>
#include <Rcpp.h>

#include <limits> /* for numeric_limits */

#include "ints.h"

/*************** TYPES      *******************/

typedef int_fast64_t cost;
typedef int16 lap_row;
typedef int16 lap_col;

#define ALL_ONES splitbit((std::numeric_limits<splitbit>::max)())

/* For a reason I've not determined, shrinking BIG is necessary to avoid 
 * an infinite loop in lap. */
#define BIG cost((std::numeric_limits<cost>::max)() / SL_MAX_SPLITS)

#define ROUND_PRECISION cost(2048 * 2048)

/***** Constants requiring initialization *****/

extern double
  lg2[int32(SL_MAX_TIPS - 1) * (SL_MAX_TIPS - 1) + 1],
  lg2_double_factorial[SL_MAX_TIPS + SL_MAX_TIPS - 2],
  lg2_unrooted[SL_MAX_TIPS + 2];
extern double *lg2_rooted;

/*************** FUNCTIONS  *******************/

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
  spi_overlap(const splitbit* a_state, const splitbit* b_state,
              const int16 n_tips, const int16 in_a, const int16 in_b,
              const int16 n_bins);

extern Rcpp::List cpp_robinson_foulds_distance(Rcpp::RawMatrix x,
                                               Rcpp::RawMatrix y,
                                               Rcpp::IntegerVector nTip);
#endif