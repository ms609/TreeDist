/************************************************************************
 *
 *  lap.cpp
 version 1.0 - 4 September 1996
 author: Roy Jonker @ MagicLogic Optimization Inc.
 e-mail: roy_jonker@magiclogic.com
 
 Code for Linear Assignment Problem, according to
 
 "A Shortest Augmenting Path Algorithm for Dense and Sparse Linear
 Assignment Problems," Computing 38, 325-340, 1987
 
 by
 
 R. Jonker and A. Volgenant, University of Amsterdam.
 
 *
 CHANGED 2016-05-13 by Yong Yang(yongyanglink@gmail.com) in column reduction 
 part according to matlab version of LAPJV algorithm (Copyright (c) 2010, Yi Cao
 All rights reserved)--
 https://www.mathworks.com/matlabcentral/fileexchange/26836-lapjv-jonker-volgenant-algorithm-for-linear-assignment-problem-v3-0:
 *
 *************************************************************************/
#include "tree_distances.h"
#include <stdio.h>
#include <iostream>
using namespace std;

/* This function is the jv shortest augmenting path algorithm to solve the 
   assignment problem */
cost lap(int dim,
         cost **assigncost,
         lap_col *rowsol,
         lap_row *colsol,
         cost *u,
         cost *v)
  
  // input:
  // dim        - problem size
  // assigncost - cost matrix
  
  // output:
  // rowsol     - column assigned to row in solution
  // colsol     - row assigned to column in solution
  // u          - dual variables, row reduction numbers
  // v          - dual variables, column reduction numbers
  
{
  bool unassignedfound;
  lap_row  i, imin, numfree = 0, prvnumfree, f, i0, k, freerow, *pred, *free;
  lap_col  j, j1, j2 = 0, endofpath, last = 0, low, up, *collist, *matches;
  /* Initializing j2 and last is unnecessary, but avoids compiler warnings */
  cost min, h, umin, usubmin, v2, *d;
  
  free = new lap_row[MAX_SPLITS];       // list of unassigned rows.
  collist = new lap_col[MAX_SPLITS];    // list of columns to be scanned in various ways.
  matches = new lap_col[MAX_SPLITS];    // counts how many times a row could be assigned.
  d = new cost[MAX_SPLITS];             // 'cost-distance' in augmenting path calculation.
  pred = new lap_row[MAX_SPLITS];       // row-predecessor of column in augmenting/alternating path.
  
  // init how many times a row will be assigned in the column reduction.
  for (i = 0; i < dim; i++)
    matches[i] = 0;
  
  // COLUMN REDUCTION
  for (j = dim; j--; ) // reverse order gives better results.
  {
    // find minimum cost over rows.
    min = assigncost[0][j];
    imin = 0;
    for (i = 1; i < dim; i++)
      if (assigncost[i][j] < min)
      {
        min = assigncost[i][j];
        imin = i;
      }
      v[j] = min;
      if (++matches[imin] == 1)
      {
        // init assignment if minimum row assigned for first time.
        rowsol[imin] = j;
        colsol[j] = imin;
      }
      else if(v[j] < v[rowsol[imin]]) {
        int j1 = rowsol[imin];
        rowsol[imin] = j;
        colsol[j] = imin;
        colsol[j1] = -1;
      }
      else
      {
        colsol[j] = -1;        // row already assigned, column not assigned.
      }
  }
  
  // REDUCTION TRANSFER
  for (i = 0; i < dim; i++)
    if (matches[i] == 0)     // fill list of unassigned 'free' rows.
      free[numfree++] = i;
    else
      if (matches[i] == 1)   // transfer reduction from rows that are assigned once.
      {
        j1 = rowsol[i];
        min = BIG;
        for (j = 0; j < dim; j++)
          if (j != j1)
            if (assigncost[i][j] - v[j] < min)
              min = assigncost[i][j] - v[j];
            v[j1] = v[j1] - min;
      }
      
      //   AUGMENTING ROW REDUCTION
      int loopcnt = 0;           // do-loop to be done twice.
      do
      {
        loopcnt++;
        
        //     scan all free rows.
        //     in some cases, a free row may be replaced with another one to be scanned next.
        k = 0;
        prvnumfree = numfree;
        numfree = 0;             // start list of rows still free after augmenting row reduction.
        while (k < prvnumfree)
        {
          i = free[k];
          k++;
          
          //       find minimum and second minimum reduced cost over columns.
          umin = assigncost[i][0] - v[0];
          j1 = 0;
          usubmin = BIG;
          for (j = 1; j < dim; j++)
          {
            h = assigncost[i][j] - v[j];
            if (h < usubmin) {
              if (h >= umin)
              {
                usubmin = h;
                j2 = j;
              }
              else
              {
                usubmin = umin;
                umin = h;
                j2 = j1;
                j1 = j;
              }
            }
          }
          
          i0 = colsol[j1];
          if (umin < usubmin) {
            //         change the reduction of the minimum column to increase the minimum
            //         reduced cost in the row to the subminimum.
            v[j1] = v[j1] - (usubmin - umin);
          }
          else {                   // minimum and subminimum equal.
            if (i0 > -1)  // minimum column j1 is assigned.
            {
              //           swap columns j1 and j2, as j2 may be unassigned.
              j1 = j2;
              i0 = colsol[j2];
            }
          }
            
          //       (re-)assign i to j1, possibly de-assigning an i0.
          rowsol[i] = j1;
          colsol[j1] = i;
          
          if (i0 > -1) { // minimum column j1 assigned earlier.
            if (umin < usubmin) {
              //           put in current k, and go back to that k.
              //           continue augmenting path i - j1 with i0.
              free[--k] = i0;
            }
            else
            {
              //           no further augmenting reduction possible.
              //           store i0 in list of free rows for next phase.
              free[numfree++] = i0;
            }
          }
        }
      }
      while (loopcnt < 2);       // repeat once.
      
      // AUGMENT SOLUTION for each free row.
      for (f = 0; f < numfree; f++)
      {
        freerow = free[f];       // start row of augmenting path.
        
        // Dijkstra shortest path algorithm.
        // runs until unassigned column added to shortest path tree.
        for(j = dim;j--;)
        {
          d[j] = assigncost[freerow][j] - v[j];
          pred[j] = freerow;
          collist[j] = j;        // init column list.
        }
        
        low = 0; // columns in 0..low-1 are ready, now none.
        up = 0;  // columns in low..up-1 are to be scanned for current minimum, now none.
        // columns in up..dim-1 are to be considered later to find new minimum,
        // at this stage the list simply contains all columns
        unassignedfound = false;
        do
        {
          if (up == low)         // no more columns to be scanned for current minimum.
          {
            last = low - 1;
            
            // scan columns for up..dim-1 to find all indices for which new minimum occurs.
            // store these indices between low..up-1 (increasing up).
            min = d[collist[up++]];
            for (k = up; k < dim; k++)
            {
              j = collist[k];
              h = d[j];
              if (h <= min)
              {
                if (h < min)     // new minimum.
                {
                  up = low;      // restart list at index low.
                  min = h;
                }
                // new index with same minimum, put on undex up, and extend list.
                collist[k] = collist[up];
                collist[up++] = j;
              }
            }
            // check if any of the minimum columns happens to be unassigned.
            // if so, we have an augmenting path right away.
            for (k = low; k < up; k++)
              if (colsol[collist[k]] < 0)
              {
                endofpath = collist[k];
                unassignedfound = true;
                break;
              }
          }
          
          if (!unassignedfound)
          {
            // update 'distances' between freerow and all unscanned columns,
            // via next scanned column.
            j1 = collist[low];
            low++;
            i = colsol[j1];
            h = assigncost[i][j1] - v[j1] - min;
            
            for (k = up; k < dim; k++)
            {
              j = collist[k];
              v2 = assigncost[i][j] - v[j] - h;
              if (v2 < d[j])
              {
                pred[j] = i;
                if (v2 == min)   // new column found at same minimum value
                {
                  if (colsol[j] < 0)
                  {
                    // if unassigned, shortest augmenting path is complete.
                    endofpath = j;
                    unassignedfound = true;
                    break;
                  }
                  // else add to list to be scanned right away.
                  else
                  {
                    collist[k] = collist[up];
                    collist[up++] = j;
                  }
                }
                d[j] = v2; // <MS: Unintended>
              }
            }
          }
        }
        while (!unassignedfound);
        
        // update column prices.
        for( k = last+1;k--;)
        {
          j1 = collist[k];
          v[j1] = v[j1] + d[j1] - min;
        }
        
        // reset row and column assignments along the alternating path.
        do
        {
          i = pred[endofpath];
          colsol[endofpath] = i;
          j1 = endofpath;
          endofpath = rowsol[i];
          rowsol[i] = j1;
        }
        while (i != freerow);
      }
      
      // calculate optimal cost.
      cost lapcost = 0;
      //  for (i = 0; i < dim; i++)
      for(i = dim;i--;)
      {
        j = rowsol[i];
        u[i] = assigncost[i][j] - v[j];
        lapcost = lapcost + assigncost[i][j];
      }
      
      // free reserved memory.
      delete[] pred;
      delete[] free;
      delete[] collist;
      delete[] matches;
      delete[] d;
      return lapcost;
}
