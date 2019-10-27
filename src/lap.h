#include "tree_distances.h" /* for score_t */

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
  typedef score_t cost;

/*************** FUNCTIONS  *******************/

extern int lap(int dim, int **assigncost,
               int *rowsol, int *colsol, int *u, int *v);

extern void checklap(int dim, int **assigncost,
                     int *rowsol, int *colsol, int *u, int *v);

