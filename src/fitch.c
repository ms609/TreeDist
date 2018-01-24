//
//  fitch.c
//  MorPhy2
//
//  Created by mbrazeau on 02/05/2017.
//  Copyright © 2017 brazeaulab. All rights reserved.
//
#include "mpl.h"
#include "morphydefs.h"
#include "morphy.h"
#include "mplerror.h"
#include "fitch.h"
#include "statedata.h"

/**/
int mpl_fitch_downpass
(MPLndsets* lset, MPLndsets* rset, MPLndsets* nset, MPLpartition* part)
{
    int i     = 0;
    int j     = 0;
    int steps = 0;
    const int* indices  = part->charindices;
    int nchars          = part->ncharsinpart;
    MPLstate* left      = lset->downpass1;
    MPLstate* right     = rset->downpass1;
    MPLstate* n         = nset->downpass1;
    
    unsigned long* weights = part->intwts;
    
    for (i = 0; i < nchars; ++i) {
        j = indices[i];
        
        n[j] = left[j] & right[j];
        
        if (n[j] == 0) {
            n[j] = left[j] | right[j];
            steps += weights[i];
        }
    }
    
    // TODO: rewrite for updated stateset checks.
    
    return steps;
}


int mpl_fitch_uppass
(MPLndsets* lset, MPLndsets* rset, MPLndsets* nset, MPLndsets* ancset,
 MPLpartition* part)
{
    int i     = 0;
    int j     = 0;
    const int* indices  = part->charindices;
    int nchars      = part->ncharsinpart;
    MPLstate* left  = lset->downpass1;
    MPLstate* right = rset->downpass1;
    MPLstate* npre  = nset->downpass1;
    MPLstate* nfin  = nset->uppass1;
    MPLstate* anc   = ancset->uppass1;
    
    for (i = 0; i < nchars; ++i) {
        
        j = indices[i];
        
        nfin[j] = anc[j] & npre[j];
        
        if (nfin[j] != anc[j]) {
            
            if (left[j] & right[j]) {
                nfin[j] = (npre[j] | (anc[j] & (left[j] | right[j])));
            }
            else {
                nfin[j] = npre[j] | anc[j];
            }
        }
#ifdef DEBUG
        assert(nfin[j]);
#endif
    }
    
    return 0;
}


int mpl_fitch_local_reopt
(MPLndsets* srcset, MPLndsets* tgt1set, MPLndsets* tgt2set, MPLpartition* part,
 int maxlen, bool domaxlen)
{
   
    int i     = 0;
    int j     = 0;
    int steps = 0;
    const int* indices    = part->charindices;
    int nchars      = part->ncharsinpart;
    MPLstate* tgt1  = tgt1set->uppass1;
    MPLstate* tgt2  = tgt2set->uppass1;
    MPLstate* src   = srcset->downpass1;
    
    unsigned long* weights = part->intwts;
    
    for (i = 0; i < nchars; ++i) {
        
        j = indices[i];
        
        if (!(src[j] & (tgt1[j] | tgt2[j]))) {
            
            steps += weights[i];
            
            if (steps > maxlen && domaxlen == true)
            {
                return steps;
            }
        }
    }
    
    return steps;
}


int mpl_NA_fitch_first_downpass
(MPLndsets* lset, MPLndsets* rset, MPLndsets* nset, MPLpartition* part)
{
    int i               = 0;
    int j               = 0;
    const int* indices  = part->charindices;
    int nchars          = part->ncharsinpart;
    MPLstate* left      = lset->downpass1;
    MPLstate* right     = rset->downpass1;
    MPLstate* n         = nset->downpass1;
    MPLstate* nt        = nset->temp_downpass1;
    
    for (i = 0; i < nchars; ++i) {
        
        j = indices[i];
        
        nset->changes[j] = false;
        
        n[j] = (left[j] & right[j]);
        
        if (n[j] == 0) {
            n[j] = (left[j] | right[j]);
            
            if ((left[j] & ISAPPLIC) && (right[j] & ISAPPLIC)) {
                n[j] = n[j] & ISAPPLIC;
            }
        }
        else {
            if (n[j] == NA) {
                if ((left[j] & ISAPPLIC) && (right[j] & ISAPPLIC)) {
                    n[j] = (left[j] | right[j]);
                }
            }
        }
        
        nt[j] = n[j]; // Store a copy for partially reoptimising the subtree
#ifdef DEBUG
        assert(n[j]);
#endif
    }
    
    return 0;
}

int mpl_nadown1_simpl
(MPLndsets* lset, MPLndsets* rset, MPLndsets* nset, MPLpartition* part)
{
    int i               = 0;
    int j               = 0;
    const int* indices  = part->charindices;
    int nchars          = part->ncharsinpart;
    MPLstate* left      = lset->downpass1;
    MPLstate* right     = rset->downpass1;
    MPLstate* n         = nset->downpass1;
    MPLstate* nt        = nset->temp_downpass1;
    
    for (i = 0; i < nchars; ++i) {
        
        j = indices[i];
        
        if (left[j] & ISAPPLIC && right[j] & ISAPPLIC) {
            n[j] = (left[j] | right[j]) & ISAPPLIC;
        }
        else {
            n[j] = (left[j] & right[j]);
            if (n[j] != NA) {
                n[j] = (left[j] | right[j]);
            }
        }
        
        nt[j] = n[j]; // Store a copy for partially reoptimising the subtree
#ifdef DEBUG
        assert(n[j]);
#endif
    }
    
    return 0;
}


int mpl_NA_fitch_first_update_downpass
(MPLndsets* lset, MPLndsets* rset, MPLndsets* nset, MPLpartition* part)
{
    /*------------------------------------------------------------------------*
     |  This function is for doing a partial downpass when proposing a        |
     |  subtree reinsertion during branchswapping. Its purpose is to          |
     |  (partially) correct any character state sets that are affected by     |
     |  the proposed reinsertion. It is nearly identical to its original-     |
     |  pass counterpart except that it does not overwrite the temp state     |
     |  storage.                                                              |
     *------------------------------------------------------------------------*/
    int i               = 0;
    int j               = 0;
    const int* indices  = part->update_NA_indices;
    int nchars          = part->nNAtoupdate;
    MPLstate* left      = lset->downpass1;
    MPLstate* right     = rset->downpass1;
    MPLstate* n         = nset->downpass1;
    MPLstate* ntemp     = nset->temp_downpass1;
    
    for (i = 0; i < nchars; ++i) {
        
        j = indices[i];
        
        n[j] = (left[j] & right[j]);
        
        if (n[j] == 0) {
            n[j] = (left[j] | right[j]);
            
            if ((left[j] & ISAPPLIC) && (right[j] & ISAPPLIC)) {
                n[j] = n[j] & ISAPPLIC;
            }
        }
        else {
            if (n[j] == NA) {
                if ((left[j] & ISAPPLIC) && (right[j] & ISAPPLIC)) {
                    n[j] = (left[j] | right[j]);
                }
            }
        }
        
        if (n[j] != ntemp[j]) {
            nset->updated = true;
        }
        
#ifdef DEBUG
        assert(n[j]);
#endif
    }
    
    return 0;
}


int mpl_NA_fitch_first_uppass
(MPLndsets* lset, MPLndsets* rset, MPLndsets* nset, MPLndsets* ancset,
 MPLpartition* part)
{
    int         i       = 0;
    int         j       = 0;
    const int*  indices = part->charindices;
    int         nchars  = part->ncharsinpart;
    MPLstate*   left    = lset->downpass1;
    MPLstate*   right   = rset->downpass1;
    MPLstate*   npre    = nset->downpass1;
    MPLstate*   nifin   = nset->uppass1;
    MPLstate*   anc     = ancset->uppass1;
    MPLstate*   nfint   = nset->temp_uppass1;
    
    for (i = 0; i < nchars; ++i) {
        
        j = indices[i];
        
        if (npre[j] & NA) {
            if (npre[j] & ISAPPLIC) {
                if (anc[j] == NA) {
                    nifin[j] = NA;
                }
                else {
                    nifin[j] = npre[j] & ISAPPLIC;
                }
            }
            else {
                if (anc[j] == NA) {
                    nifin[j] = NA;
                }
                else {
                    if ((left[j] | right[j]) & ISAPPLIC) {
                        nifin[j] = ((left[j] | right[j]) & ISAPPLIC);
                    }
                    else {
                        nifin[j] = NA;
                    }
                }
            }
        }
        else {
            nifin[j] = npre[j];
        }
        
        // Store the set for restoration during tree searches.
        nfint[j] = nifin[j];
        
#ifdef DEBUG
        assert(nifin[j]);
#endif
    }
    
    return 0;
}

int mpl_naupp1_simpl
(MPLndsets* lset, MPLndsets* rset, MPLndsets* nset, MPLndsets* ancset,
 MPLpartition* part)
{
    int         i       = 0;
    int         j       = 0;
    const int*  indices = part->charindices;
    int         nchars  = part->ncharsinpart;
    MPLstate*   left    = lset->downpass1;
    MPLstate*   right   = rset->downpass1;
    MPLstate*   npre    = nset->downpass1;
    MPLstate*   nifin   = nset->uppass1;
    MPLstate*   anc     = ancset->uppass1;
    MPLstate*   nfint   = nset->temp_uppass1;
    
    for (i = 0; i < nchars; ++i) {
        
        j = indices[i];
        
        // TODO: Rewrite this
        if (anc[j] == NA) {
            nifin[j] = npre[j];
            if (left[j] & right[j] & NA) {
                nifin[j] = NA;
            }
        }
        else {
            nifin[j] = npre[j] & anc[j];
            
            if (nifin[j] != anc[j]) {
                nifin[j] = (npre[j] | (anc[j] & (left[j] | right[j])));
            }
            else {
                nifin[j] = nifin[j] | anc[j];
            }
        }

        // Store the set for restoration during tree searches.
        nfint[j] = nifin[j];
        
#ifdef DEBUG
        assert(nifin[j]);
#endif
    }
    
    return 0;
}


int mpl_NA_fitch_first_update_uppass
(MPLndsets* lset, MPLndsets* rset, MPLndsets* nset, MPLndsets* ancset,
 MPLpartition* part)
{
    /*------------------------------------------------------------------------*
     |  This function is for doing a partial uppsass when proposing a subtree |
     |  reinsertion during branchswapping. Its purpose is to (partially)      |
     |  correct any character state sets that are affected by the proposed    |
     |  reinsertion.                                                          |
     *------------------------------------------------------------------------*/
    int         i       = 0;
    int         j       = 0;
    const int*  indices = part->update_NA_indices;
    int         nchars  = part->nNAtoupdate;
    MPLstate*   left    = lset->downpass1;
    MPLstate*   right   = rset->downpass1;
    MPLstate*   npre    = nset->downpass1;
    MPLstate*   nifin   = nset->uppass1;
    MPLstate*   anc     = ancset->uppass1;
    
    for (i = 0; i < nchars; ++i) {
        
        j = indices[i];
        
        if (npre[j] & NA) {
            if (npre[j] & ISAPPLIC) {
                if (anc[j] == NA) {
                    nifin[j] = NA;
                }
                else {
                    nifin[j] = npre[j] & ISAPPLIC;
                }
            }
            else {
                if (anc[j] == NA) {
                    nifin[j] = NA;
                }
                else {
                    if ((left[j] | right[j]) & ISAPPLIC) {
                        nifin[j] = ((left[j] | right[j]) & ISAPPLIC);
                    }
                    else {
                        nifin[j] = NA;
                    }
                }
            }
        }
        else {
            nifin[j] = npre[j];
        }
        
#ifdef DEBUG
        assert(nifin[j]);
#endif
    }
    
    return 0;
}


int mpl_NA_fitch_second_downpass
(MPLndsets* lset, MPLndsets* rset, MPLndsets* nset, MPLpartition* part)
{
    int             i       = 0;
    int             j       = 0;
    int             steps   = 0;
    const int*      indices = part->charindices;
    int             nchars  = part->ncharsinpart;
    MPLstate*       left    = lset->downpass2;
    MPLstate*       right   = rset->downpass2;
    MPLstate*       nifin   = nset->uppass1;
    MPLstate*       npre    = nset->downpass2;
    MPLstate*       npret   = nset->temp_downpass2;
    MPLstate*       stacts  = nset->subtree_actives;
    MPLstate*       tstatcs = nset->temp_subtr_actives;
    MPLstate*       lacts   = lset->subtree_actives;
    MPLstate*       racts   = rset->subtree_actives;
    MPLstate        temp    = 0;
    unsigned long*  weights = part->intwts;
    
    for (i = 0; i < nchars; ++i) {
        
        j = indices[i];
        
        nset->changes[j] = false;
        
        if (nifin[j] & ISAPPLIC) {
            if ((temp = (left[j] & right[j]))) {
                if (temp & ISAPPLIC) {
                    npre[j] = temp & ISAPPLIC;
                } else {
                    npre[j] = temp;
                }
            }
            else {
                npre[j] = (left[j] | right[j]) & ISAPPLIC;
                
                if (left[j] & ISAPPLIC && right[j] & ISAPPLIC) {
                    steps += weights[i];
                    nset->changes[j] = true;
                } else if (lacts[j] && racts[j]) {
                    steps += weights[i];
                    nset->changes[j] = true;
                }
            }
        }
        else {
            npre[j] = nifin[j];
            
            if (lacts[j] && racts[j]) {
                steps += weights[i];
                nset->changes[j] = true;
            }
        }
        
        /* Store the states active on this subtree */
        stacts[j]   = (lacts[j] | racts[j]) & ISAPPLIC;
        
        npret[j]    = npre[j]; // Storage for temporary updates.
        tstatcs[j]  = stacts[j]; // Storage for temporary updates.
    
#ifdef DEBUG
        assert(npre[j]);
#endif
    }
    
    return steps;
}


int mpl_NA_fitch_second_update_downpass
(MPLndsets* lset, MPLndsets* rset, MPLndsets* nset, MPLpartition* part)
{
    /*------------------------------------------------------------------------*
     |  This function is for doing a partial downpass when proposing a        |
     |  subtree reinsertion during branchswapping. Its purpose is to          |
     |  (partially) correct any character state sets that are affected by     |
     |  the proposed reinsertion. It is nearly identical to its original-     |
     |  pass counterpart except that it does not overwrite the temp state     |
     |  storage.                                                              |
     *------------------------------------------------------------------------*/
    int             i           = 0;
    int             j           = 0;
    int             steps       = 0;
    const int*      indices     = part->update_NA_indices;
    int             nchars      = part->nNAtoupdate;
    MPLstate*       left        = lset->downpass2;
    MPLstate*       right       = rset->downpass2;
    MPLstate*       nifin       = nset->uppass1;
    MPLstate*       npre        = nset->downpass2;
    const MPLstate* npret       = nset->temp_downpass2;
    MPLstate*       stacts      = nset->subtree_actives;
    MPLstate*       lacts       = lset->subtree_actives;
    MPLstate*       racts       = rset->subtree_actives;
    MPLstate        temp        = 0;
    unsigned long*  weights     = part->intwts;
    
    for (i = 0; i < nchars; ++i) {
        
        j = indices[i];
        
        if (nifin[j] & ISAPPLIC) {
            if ((temp = (left[j] & right[j]))) {
                if (temp & ISAPPLIC) {
                    npre[j] = temp & ISAPPLIC;
                } else {
                    npre[j] = temp;
                }
            }
            else {
                npre[j] = (left[j] | right[j]) & ISAPPLIC;
                
                if (left[j] & ISAPPLIC && right[j] & ISAPPLIC) {
                    steps += weights[i];
                } else if (lacts[j] && racts[j]) {
                    steps += weights[i];
                }
            }
        }
        else {
            npre[j] = nifin[j];
        }
        
        stacts[j] = (lacts[j] | racts[j]) & ISAPPLIC;
        
        /* Flag as updated if current set is different from previous */
        if (npre[j] != npret[j]) {
            nset->updated = true;
        }
        
#ifdef DEBUG
        assert(npre[j]);
#endif
    }
    
    return steps;
}


int mpl_NA_fitch_second_uppass
(MPLndsets* lset, MPLndsets* rset, MPLndsets* nset, MPLndsets* ancset,
 MPLpartition* part)
{
    int             i       = 0;
    int             j       = 0;
    int             steps   = 0;
    const int*      indices = part->charindices;
    int             nchars  = part->ncharsinpart;
    MPLstate*       left    = lset->downpass2;
    MPLstate*       right   = rset->downpass2;
    MPLstate*       npre    = nset->downpass2;
    MPLstate*       nfin    = nset->uppass2;
    MPLstate*       nfint   = nset->temp_uppass2;
    MPLstate*       anc     = ancset->uppass2;
//    MPLstate*       lacts   = lset->subtree_actives;
//    MPLstate*       racts   = rset->subtree_actives;
//    unsigned long*  weights = part->intwts;
    
    for (i = 0; i < nchars; ++i) {
        
        j = indices[i];
        
        if (npre[j] & ISAPPLIC) {
            if (anc[j] & ISAPPLIC) {
                if ((anc[j] & npre[j]) == anc[j]) {
                    nfin[j] = anc[j] & npre[j];
                } else {
                    if (left[j] & right[j]) {
                        nfin[j] = (npre[j] | (anc[j] & (left[j] | right[j])));
                    }
                    else {
                        if ((left[j] | right[j]) & NA) {
                            if ((left[j] | right[j]) & anc[j]) {
                                nfin[j] = anc[j];
                            } else {
                                nfin[j] = (left[j] | right[j] | anc[j]) & ISAPPLIC;
                            }
                        } else {
                            nfin[j] = npre[j] | anc[j];
                        }
                    }
                }
            }
            else {
                nfin[j] = npre[j];
            }
        }
        else {
            nfin[j] = npre[j];
            
            /*if (lacts[j] && racts[j]) {
                steps += weights[i];
                nset->changes[j] = true;
            }*/
        }
        
        nfint[j] = nfin[j]; // Storage of states for undoing temp updates
#ifdef DEBUG
        assert(nfin[j]);
#endif
    }
    
    return steps;
}


int mpl_NA_fitch_second_update_uppass
(MPLndsets* lset, MPLndsets* rset, MPLndsets* nset, MPLndsets* ancset,
 MPLpartition* part)
{
    /*------------------------------------------------------------------------*
     |  This function is for doing a partial uppsass when proposing a subtree |
     |  reinsertion during branchswapping. Its purpose is to (partially)      |
     |  correct any character state sets that are affected by the proposed    |
     |  reinsertion.                                                          |
     *------------------------------------------------------------------------*/
    int             i           = 0;
    int             j           = 0;
    int             steps       = 0;
    int             step_recall = 0;
    const int*      indices     = part->update_NA_indices;
    int             nchars      = part->nNAtoupdate;
    MPLstate*       left        = lset->downpass2;
    MPLstate*       right       = rset->downpass2;
    MPLstate*       npre        = nset->downpass2;
    MPLstate*       nfin        = nset->uppass2;
    MPLstate*       nfint       = nset->temp_uppass2;
    MPLstate*       anc         = ancset->uppass2;
    MPLstate*       lacts       = lset->subtree_actives;
    MPLstate*       racts       = rset->subtree_actives;
    unsigned long*  weights     = part->intwts;
    
    for (i = 0; i < nchars; ++i) {
        
        j = indices[i];
        
        if (npre[j] & ISAPPLIC) {
            if (anc[j] & ISAPPLIC) {
                if ((anc[j] & npre[j]) == anc[j]) {
                    nfin[j] = anc[j] & npre[j];
                } else {
                    if (left[j] & right[j]) {
                        nfin[j] = (npre[j] | (anc[j] & (left[j] | right[j])));
                    }
                    else {
                        if ((left[j] | right[j]) & NA) {
                            if ((left[j] | right[j]) & anc[j]) {
                                nfin[j] = anc[j];
                            } else {
                                nfin[j] = (left[j] | right[j] | anc[j]) & ISAPPLIC;
                            }
                        } else {
                            nfin[j] = npre[j] | anc[j];
                        }
                    }
                }
            }
            else {
                nfin[j] = npre[j];
            }
        }
        else {
            nfin[j] = npre[j];
            
            if (lacts[j] && racts[j]) {
                steps += weights[i];
            }
        }
        
        if (nfint[j] != nfin[j]) {
            nset->updated = true;
        }
        
        if (nset->changes[j] == true) {
            step_recall += weights[i];
        }
        
        
#ifdef DEBUG
        assert(nfin[j]);
#endif
    }
    
    nset->steps_to_recall += step_recall;
    
    return steps;
}


int mpl_fitch_NA_local_reopt
(MPLndsets* srcset, MPLndsets* tgt1set, MPLndsets* tgt2set, MPLpartition* part,
 int maxlen, bool domaxlen)
{
    
    part->ntoupdate = 0; // V. important: resets the record of characters needing updates
    
    int i           = 0;
    int j           = 0;
    int need_update = 0;
    int steps       = 0;
    const int* indices  = part->charindices;
    int nchars          = part->ncharsinpart;
   
    for (i = 0; i < nchars; ++i) {
        
        j = indices[i];
        
        part->update_NA_indices[need_update] = j;
        ++need_update;
        
//        if (!(src[j] & (tgt1f[j] | tgt2f[j]))) {
//            
//            if (src[j] & ISAPPLIC) {
//                if ((tgt1f[j] | tgt2f[j]) & ISAPPLIC) {
//                    steps += weights[i];
//                }
//                else {
//                    
//                    /* NOTE: This will be written simply at first, but there are
//                     * possible additional checks on tgt preliminary sets that 
//                     * could reduce the number of characters that need to be 
//                     * updated. */
//                    
//                    part->update_NA_indices[need_update] = j;
//                    ++need_update;
//                }
//            }
//            else {
////                This is a conservative omission for now.
////                if (src[j] & (tgt1d1[j] | tgt2d1[j])) {
//                    part->update_NA_indices[need_update] = j;
//                    ++need_update;
////                }
//            }
//        }
    }
    
    part->nNAtoupdate = need_update;
    
    return steps;
}


int mpl_fitch_tip_update(MPLndsets* tset, MPLndsets* ancset, MPLpartition* part)
{
    int i     = 0;
    int j     = 0;
    int* indices    = part->charindices;
    int nchars      = part->ncharsinpart;
    // TODO: Check these!!!!!!!
    MPLstate* tprelim = tset->downpass1;
    MPLstate* tfinal  = tset->uppass1;
    MPLstate* ttfinal = tset->temp_uppass1;
    MPLstate* astates = ancset->uppass1;
    
    for (i = 0; i < nchars; ++i) {
        j = indices[i];
        if (tprelim[j] & astates[j]) {
            tfinal[j] = tprelim[j] & astates[j];
        }
        else {
            tfinal[j] = tprelim[j];
        }
        ttfinal[j] = tfinal[j];
#ifdef DEBUG
        assert(tfinal[j]);
#endif
    }
    return 0;
}

int mpl_fitch_one_branch
(MPLndsets* tipanc, MPLndsets* node, MPLpartition* part)
{
    int i     = 0;
    int j     = 0;
    int* indices     = part->charindices;
    int nchars       = part->ncharsinpart;
    MPLstate* tipset = tipanc->downpass1;
    MPLstate* tipfin = tipanc->uppass1;
    MPLstate* ndset  = node->downpass1;
    MPLstate temp    = 0;
    unsigned long* weights = part->intwts;
    int length = 0;
    
    for (i = 0; i < nchars; ++i) {
        j = indices[i];
        
        temp = tipset[j] & ndset[j];

        if (temp == 0) {
            tipfin[j] = tipset[j];
            length += weights[i];
            node->uppass1[j] = ndset[j];
        }
        else {
            tipfin[j] = temp;
            node->uppass1[j] = temp;
        }
    }
    
    return length;
}


int mpl_fitch_NA_first_one_branch
(MPLndsets* tipanc, MPLndsets* node, MPLpartition* part)
{
    int i     = 0;
    int j     = 0;
    int* indices     = part->charindices;
    int nchars       = part->ncharsinpart;
    MPLstate* tipset = tipanc->downpass1;
    MPLstate* tipifin = tipanc->uppass1;
    MPLstate* ndset  = node->downpass1;
    MPLstate  temp   = 0;
    
    for (i = 0; i < nchars; ++i) {
        
        j = indices[i];
        
        tipanc->changes[j] = false;
        temp = tipset[j] & ndset[j];
        
        if (temp != 0) {
            tipifin[j]          = temp;
            node->uppass1[j]    = temp;
        }
    }
    
    return 0;
}

int mpl_fitch_NA_second_one_branch
(MPLndsets* tipanc, MPLndsets* node, MPLpartition* part)
{
    int i     = 0;
    int j     = 0;
    int* indices            = part->charindices;
    int nchars              = part->ncharsinpart;
    MPLstate* tipset        = tipanc->downpass1;
    MPLstate* tipifin       = tipanc->uppass1;
    MPLstate* ndset         = node->downpass2;
    MPLstate* ndacts        = node->subtree_actives;
    MPLstate  temp          = 0;
    unsigned long* weights  = part->intwts;
    int length = 0;
    
    for (i = 0; i < nchars; ++i) {
        
        j = indices[i];
        
        temp = tipset[j] & ndset[j];
        
        if (temp == 0) {
            if (tipset[j] & ISAPPLIC) {
                if (ndset[j] & ISAPPLIC) {
                    length += weights[i];
                    tipanc->changes[j] = true;
                }
                else {
                    if (ndacts[j]) {
                        length += weights[i];
                        tipanc->changes[j] = true;
                    }
                }
            }
            
            tipifin[j]        = tipset[j];
        }
        else {
            tipifin[j]        = temp;
        }
        
        tipanc->temp_downpass1[j]   = tipanc->downpass1[j];
        tipanc->temp_uppass1[j]     = tipanc->uppass1[j];
        tipanc->temp_downpass2[j]   = tipanc->downpass2[j];
        tipanc->temp_uppass2[j]     = tipanc->uppass2[j];
    }
    
    return length;
}

int mpl_fitch_NA_second_one_branch_recalc
(MPLndsets* tipanc, MPLndsets* node, MPLpartition* part)
{
    int i     = 0;
    int j     = 0;
    int* indices            = part->charindices;
    int nchars              = part->ncharsinpart;
    MPLstate* tipset        = tipanc->downpass1;
    MPLstate* tipifin       = tipanc->uppass1;
    MPLstate* ndset         = node->downpass2;
    MPLstate* ndacts        = node->subtree_actives;
    MPLstate  temp          = 0;
    unsigned long* weights  = part->intwts;
    unsigned long step_recall = 0;
    int length = 0;
    
    for (i = 0; i < nchars; ++i) {
        
        j = indices[i];
        
        temp = tipset[j] & ndset[j];
        
        if (temp == 0) {
            if (tipset[j] & ISAPPLIC) {
                if (ndset[j] & ISAPPLIC) {
                    length += weights[i];
                }
                else {
                    if (ndacts[j]) {
                        length += weights[i];
                    }
                }
            }
            
            tipifin[j]        = tipset[j];
        }
        else {
            tipifin[j]        = temp;
        }
        
        if (tipanc->changes[j] == true) {
            step_recall += weights[i];
        }
    }
    
    tipanc->steps_to_recall += step_recall;
    
    return length;
}


int mpl_fitch_NA_tip_update
(MPLndsets* tset, MPLndsets* ancset, MPLpartition* part)
{
    int i     = 0;
    int j     = 0;
    int* indices        = part->charindices;
    int nchars          = part->ncharsinpart;
    
    MPLstate* tpass1    = tset->downpass1;
    MPLstate* tpass2    = tset->uppass1;
    MPLstate* tpass3    = tset->downpass2;
    MPLstate* ttpass1   = tset->temp_downpass1;
    MPLstate* ttpass2   = tset->temp_uppass1;
    MPLstate* ttpass3   = tset->temp_downpass2;
    MPLstate* astates   = ancset->uppass1;
    MPLstate* stacts    = tset->subtree_actives;
    MPLstate* tstatcs   = tset->temp_subtr_actives;
    
    for (i = 0; i < nchars; ++i) {
        
        j = indices[i];
        
        if (tpass1[j] & astates[j]) {
            stacts[j] = (tpass1[j] & astates[j] & ISAPPLIC);
        }
        else {
            stacts[j] |= tpass1[j] & ISAPPLIC;
        }

        tpass2[j] = tpass1[j];
        
        if (tpass2[j] & astates[j]) {
            if (astates[j] & ISAPPLIC) {
                tpass2[j] &= ISAPPLIC;
            }
        }
        
        tpass3[j]  = tpass2[j];
        
        // Store the temp sets for restoring after temporary updates
        ttpass1[j] = tpass1[j];
        ttpass2[j] = tpass2[j];
        ttpass3[j] = tpass3[j];
        tstatcs[j] = stacts[j];
#ifdef DEBUG   
        assert(tpass3[j]);
        assert(tpass2[j]);
#endif
    }
    
    return 0;
}


int mpl_fitch_NA_tip_recalc_update
(MPLndsets* tset, MPLndsets* ancset, MPLpartition* part)
{
    int i               = 0;
    int j               = 0;
    int* indices        = part->charindices;
    int nchars          = part->ncharsinpart;
    
    MPLstate* tpass1    = tset->downpass1;
    MPLstate* tpass2    = tset->uppass1;
    MPLstate* tpass3    = tset->downpass2;
    MPLstate* astates   = ancset->uppass1;
    MPLstate* stacts    = tset->subtree_actives;
    
    for (i = 0; i < nchars; ++i) {
        
        j = indices[i];
        
        if (tpass1[j] & astates[j]) {
            stacts[j] = (tpass1[j] & astates[j] & ISAPPLIC);
        }
        else {
            stacts[j] |= tpass1[j] & ISAPPLIC;
        }

        tpass2[j] = tpass1[j];
        
        if (tpass2[j] & astates[j]) {
            if (astates[j] & ISAPPLIC) {
                tpass2[j] &= ISAPPLIC;
            }
        }
        
        tpass3[j]  = tpass2[j];
        
#ifdef DEBUG   
        assert(tpass3[j]);
        assert(tpass2[j]);
#endif
    }
    
    return 0;
}


int mpl_fitch_NA_tip_finalize
(MPLndsets* tset, MPLndsets* ancset, MPLpartition* part)
{
    int i     = 0;
    int j     = 0;
    int* indices    = part->charindices;
    int nchars      = part->ncharsinpart;
    MPLstate* tpass1    = tset->downpass1;
    MPLstate* tfinal    = tset->uppass2;
    MPLstate* ttfinal   = tset->temp_uppass2;
    MPLstate* astates   = ancset->uppass2;
    MPLstate* stacts    = tset->subtree_actives;
    MPLstate* tstacts   = tset->temp_subtr_actives;
    
    for (i = 0; i < nchars; ++i) {
        
        j = indices[i];
        
        if (tpass1[j] & astates[j]) {
            tfinal[j] = tpass1[j] & astates[j];
        }
        else {
            tfinal[j] = tpass1[j];
        }
        
        //stacts[j] = tfinal[j] & ISAPPLIC;
        
        // Store the temp buffers:
        ttfinal[j] = tfinal[j];
        tstacts[j] = stacts[j];
#ifdef DEBUG
        assert(tfinal[j]);
#endif
    }
    
    return 0;
}
