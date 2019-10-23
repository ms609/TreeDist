/************************************************************************
*
*  lap.h
   version 1.0 - 21 june 1996
   author  Roy Jonker, MagicLogic Optimization Inc.
   
   header file for LAP
*
**************************************************************************/

/*************** CONSTANTS  *******************/

  /* Increased by factor of 10 from original to improve precision.
   * This allows 4294 BIG scores to be added before reaching 2^32.
   * If the maximum n_tips is increased from 3200, BIG may need to shrink
   * to avoid an overflow.  */
  #define BIG 1000000

/*************** TYPES      *******************/

  typedef int lap_row;
  typedef int lap_col;
  typedef int cost;

/*************** FUNCTIONS  *******************/

extern int lap(int dim, int **assigncost,
               int *rowsol, int *colsol, int *u, int *v);

extern void checklap(int dim, int **assigncost,
                     int *rowsol, int *colsol, int *u, int *v);

