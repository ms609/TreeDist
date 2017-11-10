//
//  fitch.c
//  MorPhy2
//
//  Created by mbrazeau on 02/05/2017.
//  Copyright © 2017 brazeaulab. All rights reserved.
//
#include "morphydefs.h"
#include "mplerror.h"
#include "morphy.h"
#include "fitch.h"
#include "statedata.h"

/**/
int mpl_fitch_downpass
(MPLndsets* lset, MPLndsets* rset, MPLndsets* nset, MPLpartition* part)
{
    int i     = 0;
    int j     = 0;
    int steps = 0;
    int* indices    = part->charindices;
    int nchars      = part->ncharsinpart;
    MPLstate* left  = lset->downpass1;
    MPLstate* right = rset->downpass1;
    MPLstate* n     = nset->downpass1;
    
    unsigned long* weights = part->intwts;
    
    for (i = 0; i < nchars; ++i) {
        j = indices[i];
        if (left[j] & right[j]) {
            n[j] = left[j] & right[j];
        }
        else {
            n[j] = left[j] | right[j];
            steps += weights[i];
        }
    }
    
    return steps;
}


int mpl_fitch_uppass
(MPLndsets* lset, MPLndsets* rset, MPLndsets* nset, MPLndsets* ancset,
 MPLpartition* part)
{
    int i     = 0;
    int j     = 0;
    int* indices    = part->charindices;
    int nchars      = part->ncharsinpart;
    MPLstate* left  = lset->downpass1;
    MPLstate* right = rset->downpass1;
    MPLstate* npre  = nset->downpass1;
    MPLstate* nfin  = nset->uppass1;
    MPLstate* anc   = ancset->uppass1;
    
    for (i = 0; i < nchars; ++i) {
        
        j = indices[i];
        
        if ((anc[j] & npre[j]) == anc[j]) {
            nfin[j] = anc[j] & npre[j];
        }
        else {
            if (left[j] & right[j]) {
                nfin[j] = (npre[j] | (anc[j] & (left[j] | right[j])));
            }
            else {
                nfin[j] = npre[j] | anc[j];
            }
        }
       
        assert(nfin[j]);
    }
    
    return 0;
}


/**/
int mpl_NA_fitch_first_downpass
(MPLndsets* lset, MPLndsets* rset, MPLndsets* nset, MPLpartition* part)
{
    int i     = 0;
    int j     = 0;
    int* indices    = part->charindices;
    int nchars      = part->ncharsinpart;
    MPLstate* left  = lset->downpass1;
    MPLstate* right = rset->downpass1;
    MPLstate* n     = nset->downpass1;
    MPLstate* stacts  = nset->subtree_actives;
    MPLstate* lacts   = lset->subtree_actives;
    MPLstate* racts   = rset->subtree_actives;
    MPLstate temp = 0;
    
    for (i = 0; i < nchars; ++i) {
        j = indices[i];
        
        if ((temp = (left[j] & right[j]))) {
            n[j] = temp;
            
            if (temp == NA) {
                if ((left[j] & ISAPPLIC) && (right[j] & ISAPPLIC)) {
                    n[j] = (left[j] | right[j]);
                }
            }
        }
        else {
            n[j] = (left[j] | right[j]);
            
            if ((left[j] & ISAPPLIC) && (right[j] & ISAPPLIC)) {
                n[j] = n[j] & ISAPPLIC;
            }
        }
        
        stacts[j] = (lacts[j] | racts[j]) & ISAPPLIC;
        
        assert(n[j]);
    }
    
    return 0;
}


int mpl_NA_fitch_first_uppass
(MPLndsets* lset, MPLndsets* rset, MPLndsets* nset, MPLndsets* ancset,
 MPLpartition* part)
{
    int i     = 0;
    int j     = 0;
    int* indices    = part->charindices;
    int nchars      = part->ncharsinpart;
    MPLstate* left  = lset->downpass1;
    MPLstate* right = rset->downpass1;
    MPLstate* npre  = nset->downpass1;
    MPLstate* nifin = nset->uppass1;
    MPLstate* anc   = ancset->uppass1;
    
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
        
        assert(nifin[j]);
    }
    
    
    return 0;
}


int mpl_NA_fitch_second_downpass
(MPLndsets* lset, MPLndsets* rset, MPLndsets* nset, MPLpartition* part)
{
    int i     = 0;
    int j     = 0;
    int steps = 0;
    int* indices    = part->charindices;
    int nchars      = part->ncharsinpart;
    MPLstate* left  = lset->downpass2;
    MPLstate* right = rset->downpass2;
    MPLstate* nifin = nset->uppass1;
    MPLstate* npre    = nset->downpass2;
    MPLstate* stacts  = nset->subtree_actives;
    MPLstate* lacts   = lset->subtree_actives;
    MPLstate* racts   = rset->subtree_actives;
    MPLstate temp = 0;
    
    unsigned long* weights = part->intwts;
    
    for (i = 0; i < nchars; ++i) {
        
//        temp = 0;
        
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
    
        assert(npre[j]);
    }
    
    return steps;
}


int mpl_NA_fitch_second_uppass
(MPLndsets* lset, MPLndsets* rset, MPLndsets* nset, MPLndsets* ancset,
 MPLpartition* part)
{
    int i     = 0;
    int j     = 0;
    int steps = 0;
    int* indices    = part->charindices;
    int nchars      = part->ncharsinpart;
    MPLstate* left  = lset->downpass2;
    MPLstate* right = rset->downpass2;
    MPLstate* npre  = nset->downpass2;
    MPLstate* nfin  = nset->uppass2;
    MPLstate* anc   = ancset->uppass2;
    MPLstate* lacts = lset->subtree_actives;
    MPLstate* racts = rset->subtree_actives;
    unsigned long* weights = part->intwts;
    
    for (i = 0; i < nchars; ++i) {
        
        j = indices[i];
        
        if (npre[j] & ISAPPLIC) {
            if (anc[j] & ISAPPLIC) {
                if ((anc[j] & npre[j]) == anc[j]) {
                    nfin[j] = anc[j] & npre[j];
                } else {
                    if (left[j] & right[j]) {
                        nfin[j] = (npre[j] | (anc[j] & left[j] & right[j]));
                    }
                    else {
                        if ((left[j] | right[j]) & NA) {
                            if ((left[j] | right[j]) & anc[j]) {
                                nfin[j] = anc[j];
                            } else {
                                nfin[j] = (left[j] | right[j] | anc[j]) & NA;
                            }
                        } else {
                            nfin[j] = npre[j] | anc[j];
                            if ((anc[j] & nfin[j]) == anc[j]) {
                                nfin[j] = anc[j] & nfin[j];
                            }
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
        assert(nfin[j]);
    }
    
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
    MPLstate* astates = ancset->uppass1;
    
    for (i = 0; i < nchars; ++i) {
        j = indices[i];
        if (tprelim[j] & astates[j]) {
            tfinal[j] = tprelim[j] & astates[j];
        }
        else {
            tfinal[j] = tprelim[j];
        }
        assert(tfinal[j]);
    }
    return 0;
}

int mpl_fitch_NA_tip_update
(MPLndsets* tset, MPLndsets* ancset, MPLpartition* part)
{
    int i     = 0;
    int j     = 0;
    int* indices    = part->charindices;
    int nchars      = part->ncharsinpart;
    MPLstate* tpass1    = tset->downpass1;
    MPLstate* tpass2   = tset->uppass1;
    MPLstate* tpass3   = tset->downpass2;
//    MPLstate* tifinal   = tset->uppass2;
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
        
        tpass3[j] = tpass2[j];
        
        assert(tpass3[j]);
        assert(tpass2[j]);
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
    MPLstate* astates   = ancset->uppass2;
    MPLstate* stacts    = tset->subtree_actives;
    
    for (i = 0; i < nchars; ++i) {
        
        j = indices[i];
        
        if (tpass1[j] & astates[j]) {
            tfinal[j] = tpass1[j] & astates[j];
        }
        else {
            tfinal[j] = tpass1[j];
        }
        
        stacts[j] = tfinal[j] & ISAPPLIC;
        assert(tfinal[j]);
    }
    
    return 0;
}
