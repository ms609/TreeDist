/*
//  morphydefs.h
//  MorPhy2
//
//  Created by mbrazeau on 07/05/2017.
//  Copyright © 2017 brazeaulab. All rights reserved.
*/

#ifndef morphydefs_h
#define morphydefs_h

#ifdef __cplusplus
extern "C" {
#endif
  
#include <stdint.h>
#include <stdbool.h>
  
  
  typedef double Mflt;
#define MPL_EPSILON DBL_EPSILON
  
  typedef unsigned int MPLstate;
  
#define NA              ((MPLstate)1)
#define MISSING         ((MPLstate)~0)
#define ISAPPLIC        (((MPLstate)~0)^NA)
#define UNKNOWN         ISAPPLIC
#define MAXSTATES       (CHAR_BIT * sizeof(MPLstate))
#define DEFAULTGAP      '-'
#define DEFAULTMISSING  '?'
#define DEFAULTUNKNOWN  '+'
#define DEFAULCHARTYPE  FITCH_T
#define DEFAULTWTBASE   1
#define NACUTOFF        2
#define MPLCHARMAX      INT_MAX
#define USRWTMIN        0.00001 /*! Minimum fractional weight a caller can ask 
  for when setting weights. Anything less than 
this will be considered 0. */
#define MPLWTMIN        (MPL_EPSILON * 10) /*! Safest (for me!) if calculations
steer pretty clear of epsilon */

typedef struct MPLndsets MPLndsets;
typedef struct MPLpartition MPLpartition;

typedef int (*MPLdownfxn)
  (MPLndsets*     lset,
   MPLndsets*     rset,
   MPLndsets*     nset,
   MPLpartition*  part);

typedef int (*MPLupfxn)
  (MPLndsets*     lset,
   MPLndsets*     rset,
   MPLndsets*     nset,
   MPLndsets*     ancset,
   MPLpartition*  part);

typedef int (*MPLtipfxn)
  (MPLndsets*     tset,
   MPLndsets*     ancset,
   MPLpartition*  part);

typedef int (*MPLloclfxn)
  (MPLndsets* srcset,
   MPLndsets* topnod,
   MPLndsets* botnod,
   MPLpartition* part,
   int cutoff,
   bool usemax);


typedef struct {
  MPLstate    asint;
  char*       asstr;
} MPLcell;


typedef struct MPLcharinfo MPLcharinfo;
struct MPLcharinfo {
  
  int         charindex;
  int         ninapplics;
  
  MPLchtype   chtype;
  double      realweight;
  unsigned long        basewt;
  unsigned long        intwt;
  Mflt        fltwt;
  Mflt         CIndex;
  Mflt         RCIndex;
  Mflt         HIndex;
  Mflt         RetIndex;
  
};


typedef struct {
  int     maxchars;
  int     nupdate;
  int*    indices;
} MPLcupdate;


struct MPLpartition {
  
  MPLchtype       chtype;         /*!< The optimality type used for this partition. */
  bool            isNAtype;       /*!< This character should be treated as having inapplicable data. */ 
  int             maxnchars;
  int             ncharsinpart;
  int*            charindices;
  unsigned long   nchanges;       /*!< Number of state changes in this partition. */
  int             ntoupdate;
  int*            update_indices;
  int             nNAtoupdate;
  int*            update_NA_indices;
  bool            usingfltwt;
  unsigned long*  intwts;
  Mflt*           fltwts;
  MPLtipfxn       tipupdate;
  MPLtipfxn       tipfinalize;
  MPLtipfxn       tiproot;        /*!< For the function that adds length at the base of an unrooted tree. */
  MPLtipfxn       tiprootfinal;
  MPLtipfxn       tipupdaterecalc;
  MPLtipfxn       tipfinalrecalc;
  MPLtipfxn       tiprootrecalc;
  MPLtipfxn       tiprootupdaterecalc;
  MPLdownfxn      inappdownfxn;
  MPLdownfxn      inappdownrecalc2;
  MPLupfxn        inappupfxn;
  MPLupfxn        inapuprecalc2;
  MPLdownfxn      prelimfxn;
  MPLdownfxn      downrecalc1;
  MPLupfxn        finalfxn;
  MPLupfxn        uprecalc1;
  MPLloclfxn      loclfxn;
  MPLpartition*   next;
  
};


struct MPLndsets {
  
  bool        updated;
  int         steps_to_recall;
  MPLstate*   downpass1;
  MPLstate*   uppass1;
  MPLstate*   downpass2;
  MPLstate*   uppass2;
  MPLstate*   subtree_actives;
  MPLstate*   temp_subtr_actives;
  MPLstate*   temp_downpass1;
  MPLstate*   temp_uppass1;
  MPLstate*   temp_downpass2;
  MPLstate*   temp_uppass2;
  bool*       changes;
  char**      downp1str;
  char**      downp2str;
  char**      upp1str;
  char**      upp2str;
  
};


typedef struct mpl_matrix_s {
  int             ncells;
  MPLcell*        cells;
} MPLmatrix;


typedef struct symbols_s {
  int         numstates;
  char*       statesymbols;
  char*       symbolsinmatrix;
  MPLstate*   packed;
  char        gap;
  char        missing;
} MPLsymbols;


/*! \struct*/
typedef struct Morphy_t {
  
  int             numtaxa;    
  int             numcharacters;
  int             numrealwts;
  MPLcharinfo*    charinfo;   
  unsigned long   usrwtbase;
  unsigned long   wtbase;   
  int             numparts; 
  MPLpartition*   partstack;  
  MPLpartition**  partitions; 
  MPLsymbols      symbols;    
  MPLgap_t           gaphandl;
  union {
    int         asint;
    Mflt        asfloat;
  } score;   
  MPLmatrix       inmatrix;   
  char*           char_t_matrix;  
  int             numnodes;   
  int*            nodesequence;   
  int             nthreads; 
  MPLndsets**     statesets;
  
} Morphy_t, *Morphyp;


#ifdef __cplusplus
}
#endif

#endif /* morphydefs_h */
