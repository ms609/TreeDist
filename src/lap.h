#include <stdint.h>

/************************************************************************
*
*  lap.h
   version 1.0 - 21 june 1996
   author  Roy Jonker, MagicLogic Optimization Inc.
   
   header file for LAP
*
**************************************************************************/

/*************** CONSTANTS  *******************/

  /* BIG is now defined in tree_distances.h */

/*************** TYPES      *******************/

  typedef int lap_row;
  typedef int lap_col;
  typedef int64_t cost;

/*************** FUNCTIONS  *******************/

extern cost lap(int dim, cost **assigncost,
               lap_col *rowsol, lap_row *colsol, cost *u, cost *v);

extern void checklap(int dim, cost **assigncost,
                     lap_col *rowsol, lap_row *colsol, cost *u, cost *v);

