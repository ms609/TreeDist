//
//  morphydefs.h
//  MorPhy2
//
//  Created by mbrazeau on 07/05/2017.
//  Copyright © 2017 brazeaulab. All rights reserved.
//

#ifndef morphydefs_h
#define morphydefs_h

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>
#include <stdbool.h>
    
//#ifdef MPLDBL
typedef double Mflt;
#define MPL_EPSILON DBL_EPSILON
//#elif MPLLDBL
//typedef long double Mflt;
//#else
//typedef float Mflt;
//#endif

typedef unsigned int MPLstate;

#define NA              ((MPLstate)01) // 0b1 is not C11 compliant
#define MISSING         ((MPLstate)~0)
#define ISAPPLIC        (((MPLstate)~0)^NA)
#define MAXSTATES       (CHAR_BIT * sizeof(MPLstate))
#define DEFAULTGAP      '-'
#define DEFAULTMISSING  '?'
#define DEFAULCHARTYPE  FITCH_T
#define DEFAULTWTBASE   1
#define NACUTOFF        2   // The max number of NA tokens that can be ignored
                            // in a column
#define MPLCHARMAX      INT_MAX
#define USRWTMIN        0.00001 /*! Minimum fractional weight a caller can ask 
                                    for when setting weights. Anything less than 
                                    this will be considered 0. */
#define MPLWTMIN        (MPL_EPSILON * 10) /*! Safest (for me!) if calculations
                                               steer pretty clear of epsilon */
    
typedef struct MPLndsets MPLndsets;
typedef struct partition_s MPLpartition;
// Evaluator function pointers
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

// Key data types
typedef struct {
    MPLstate    asint;
    char*       asstr;
} MPLcell;
    

typedef struct charinfo_s MPLcharinfo;
typedef struct charinfo_s {
    
    int         charindex;
    int         ninapplics;
//    bool        included;
    MPLchtype   chtype;
    double      realweight;
    unsigned long        basewt;
    unsigned long        intwt;
    Mflt        fltwt;
//    union {
//        unsigned long   intwt;
//        Mflt            fltwt;
//    };
    Mflt         CIndex;
    Mflt         RCIndex;
    Mflt         HIndex;
    Mflt         RetIndex;
    
} MPLcharinfo;
    
    
typedef struct {
    int     maxchars;
    int     nupdate;
    int*    indices;
} MPLcupdate;
    
    
typedef struct partition_s MPLpartition;
typedef struct partition_s {
    
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
    
} MPLpartition;
    


typedef struct MPLndsets {
    
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
    
} MPLndsets;
    
    
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
    
    int             numtaxa;    // The number of terminal taxa
    int             numcharacters;  // The number of characters (transformation series)
    int             numrealwts;
    MPLcharinfo*    charinfo;   // Data type information about each character
    unsigned long   usrwtbase;
    unsigned long   wtbase;   // Used to rescale factional weights (1 by default)
    int             numparts;   // The number of data type partitions
    MPLpartition*   partstack;  // A place for unused partitions
    MPLpartition**  partitions; // The array of partitions
    MPLsymbols      symbols;    // The symbols used in the dataset
    MPLgap_t           gaphandl;   // The method of gap treatment
    union {
        int         asint;
        Mflt        asfloat;
    } score;   // The score (parsimony, likelihood etc.) of the evaluated data
    MPLmatrix       inmatrix;   // Internal representation of the matrix
    char*           char_t_matrix;  // The matrix as a NULL-terminated string
    int             numnodes;   // The number of nodes
    int*            nodesequence;   // The postorder sequence of nodes.
    int             nthreads;   // For programs that wish to multithread
    MPLndsets**     statesets;
    
} Morphy_t, *Morphyp;


#ifdef __cplusplus
}
#endif
    
#endif /* morphydefs_h */
